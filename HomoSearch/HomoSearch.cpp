#include "HomoSearch.h"
#include "vector_unordered_map.h"


HomoPairs::HomoPairs(PAIR_TYPE& base): base(base) {
    pairs = PAIR_MAP ();
}


HomoPairs::~HomoPairs() {}


void HomoPairs::addPair(const int& i, const int& j, double sim) {
    pairs[i * base + j] += sim;
}


HomoSearch::HomoSearch(const SequenceList& seqList): seqList(seqList) {}


HomoSearch::~HomoSearch() {}


LSH::LSH(const SequenceList& seqList, const GappedKmerScan& scanner, const int& mode, const double& similarityThreshold): 
    HomoSearch(seqList),
    scanner(scanner),
    mode(mode),
    pstable_hashTable_k(3),
    pstable_hashTable_L(30),
    minHash_hashTable_k(2),
    minHash_hashRable_L(50),
    similarityThreshold(similarityThreshold),
    projectionWidth(0.1),
    table_size(100000),
    location(0.0),
    scale(1.0),
    table_step(2.0 * scale / table_size) {
        std::cout << "\n=== Start homo searching based on LSH" << std::endl;

        cauchyTable = std::vector<double> (table_size + 5, 0); 
        combinations = std::vector <std::vector <int> > (11, std::vector <int> (11, 0));

        result = PAIR_MAP();
}


LSH::~LSH() {}


void LSH::makeCauchyTable() { 
    const double pi = 3.14159265358979323846;

    //Generate a look-up table for the probability density function of the Cauchy distribution
    double x = location - scale;
    for (int i = 0; i < table_size; i++) {
        double pdf = scale / (pi * (scale * scale + (x - location) * (x - location)));
        cauchyTable[i] = pdf * table_step;
        x += table_step;
    }
}

void LSH::generateCombinations() {
    for (int i = 0; i <= 10; ++i) {
        combinations[i][0] = 1;
        combinations[i][i] = 1;
    }
    for (int i = 1; i <= 10; ++i) {
        for (int j = 1; j < i; ++j) {
            combinations[i][j] = combinations[i-1][j-1] + combinations[i-1][j];
        }
    }
}


int LSH::hashFunction(const HASH_TYPE& str){
    int hash = 5381;
    int c, index = 0;

    while ((c = str[index++])) {
        hash = ((hash << 5) + hash) + c;
    }

    return hash;
}


int LSH::getFingerprint_pStable(const int& level, const std::vector <PROFILE_TYPE>& seqProf, 
    std::unordered_map <HASH_TYPE, double>& curIndice, const int& seqId, const double& projectionOffset) {
    int bucketId = 0;
    double dot_product = 0.0;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::random_device rd;
    std::mt19937 gen(rd());
    for(auto kmer: seqProf[level]) {
        if(curIndice.count(kmer.first) == 0) {
            double u = distribution(gen);
            int index = static_cast<int>(u * table_size); 
            double x = location - scale + index * table_step;
            curIndice[kmer.first] = x + cauchyTable[index] * (u - static_cast<double>(index) / table_size) * table_size;
        }
        dot_product += kmer.second * curIndice.at(kmer.first);
    }
    //std::cout << dot_product << "\t";
    int maxKmerNumber = 2 * ( combinations[scanner.windowSize[level]][scanner.kmerLength[level]] 
        + (seqList.list[seqId].sequence.size() - scanner.windowSize[level]) 
        * combinations[scanner.windowSize[level] - 1][scanner.kmerLength[level] - 1] );

    double hashValue = (dot_product / maxKmerNumber + projectionOffset) * 1.0 / projectionWidth;
    bucketId = static_cast<int>(floor(hashValue));
    return bucketId;
}


std::vector <int> LSH::getFingerprint_minHash(const int& level, const std::vector <PROFILE_TYPE>& seqProf, 
    const int& fingerprintNumber, const int& MOD) {
    std::priority_queue <int, std::vector<int>, std::greater<int>> minValue;
    std::vector <int> fingerPrint;

    for(auto kmer: seqProf[level]) {
        HASH_TYPE str = kmer.first;
        int value = hashFunction(str) % MOD;
        for(int i = 0; i < kmer.second; i++) {
            minValue.emplace(value);
        }
    }

    for(int i = 0; i < fingerprintNumber; i ++) {
        fingerPrint.emplace_back(minValue.top());
        minValue.pop();
    }
    return fingerPrint;
}


PAIR_MAP LSH::divideBuckets(const int& levelNumber, const std::vector <int>& bucket, const PAIR_MAP& pairMap) {
    if(levelNumber == scanner.levelNumber) {
        return pairMap;
    }

    PAIR_TYPE base = seqList.list.size();
    HomoPairs preliminaryPairs = HomoPairs(base);

    int hashTable_L;
    if(mode == 1) {
        hashTable_L = pstable_hashTable_L;
    } else {
        hashTable_L = minHash_hashRable_L;
    }

    for(int rep = 0; rep < hashTable_L; rep++) {
        std::vector <std::vector <int> > fingerprints(base, std::vector <int>());
        std::unordered_map <std::vector <int>, std::vector <int>, vector_hash> hashBins;

        //mode = 1, using LSH based on 1-stable; mode = 2, using minHash
        if(mode == 1) {
            for(int k = 0; k < pstable_hashTable_k; k ++) {
                std::unordered_map<HASH_TYPE, double> indices;

                auto rng = std::mt19937{std::random_device{}()};
                auto dist = std::uniform_real_distribution<>(0, projectionWidth);
                double projectionOffset = dist(rng);

                for(auto seqId: bucket) {
                    int fingerprint = getFingerprint_pStable
                        (levelNumber, scanner.profile[seqId], indices, seqId, projectionOffset);
                    fingerprints[seqId].emplace_back(fingerprint);
                }
            }
        } else {  
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<int> distribution;
            int MOD = distribution(gen);

            for(auto seqId: bucket) {
                fingerprints[seqId] = getFingerprint_minHash(levelNumber, scanner.profile[seqId], 
                    minHash_hashTable_k, MOD);
            }
        }

        for(auto seqId: bucket) {
            //std::cout << seqId << "\t";
            //for(auto f: fingerprints[seqId]) std::cout << f << "\t";
            //std::cout << std::endl;         
            hashBins[fingerprints[seqId]].emplace_back(seqId);
        }
        
        for(auto bin: hashBins) {
            int binSize = bin.second.size();
            for(int i = 0; i < binSize; i ++) {
                for(int j = i + 1; j < binSize; j ++) {
                    preliminaryPairs.addPair(bin.second[i], bin.second[j], 1.0 / hashTable_L);
                }
            }
        }
    }

    //Discrete primitive point sets
    int bucketSize = bucket.size();
    int order = 0;
    std::unordered_map <int, int> renumber;
    for(auto seqId: bucket) {
        renumber[seqId] = order++;
    }

    //Concatenate sequence pairs using union sets
    UnionFind Union = UnionFind(bucketSize);
    for(auto p: preliminaryPairs.pairs) {
        if(p.second >= similarityThreshold) {
            int renumFirst = renumber[p.first / base], renumSecond = renumber[p.first % base];
            if(!Union.connected(renumFirst, renumSecond)) {
                Union.merge(renumFirst, renumSecond);
            }
        }
    }

    //Find connected components using union sets
    std::unordered_map <int, std::vector <int> > newbBuckets;
    for(int i = 0; i < bucketSize; i ++) {
        newbBuckets[Union.find(i)].emplace_back(bucket[i]);
    }
    for(auto x: newbBuckets) {
        for(auto e: x.second) {
            std::cout << e << "\t";
        }
        std::cout << std::endl;
    }

    //Hierarchical iteration to determine if convergence is achieved
    std::cout << "+++ End executing hierarchical LSH with window size is " << scanner.windowSize[levelNumber] 
        << " and kmer length is " << scanner.kmerLength[levelNumber] << std::endl;
    std::cout << "+++ In this step " << newbBuckets.size() << " area(s) are delineated.\n" << std::endl;
    if(newbBuckets.size() == 1) {
        return preliminaryPairs.pairs;
    } else {
        PAIR_MAP totalMap;
        for(auto bucketElement: newbBuckets) {
            int bucketElementSize = bucketElement.second.size();
            if(bucketElementSize > 1) {
                PAIR_MAP subMap;
                for(int i = 0; i < bucketElementSize; i ++) {
                    for(int j = i + 1; j < bucketElementSize; j ++) {
                        PAIR_TYPE key = bucketElement.second[i] * base + bucketElement.second[j];
                        subMap[key] = preliminaryPairs.pairs[key];
                    }
                }
                PAIR_MAP tmp = divideBuckets(levelNumber + 1, bucketElement.second, subMap);
                for(auto element: tmp) {
                    totalMap.emplace(element.first, element.second);
                }
            }
        }
        return totalMap;
    }
}


void LSH::work() {
    makeCauchyTable();
    generateCombinations();

    int sequenceNumber = scanner.profile.size();
    std::vector <int> originalBucket = std::vector <int> ();
    for(int i = 0; i < sequenceNumber; i ++) {
        originalBucket.emplace_back(i);
    }
    result = divideBuckets(0, originalBucket, PAIR_MAP());
}


void LSH::dumpToFile(const std::string& outputFilename) {
    std::ofstream ofs(outputFilename);
    long long pairNum = 0;
    PAIR_TYPE base = seqList.list.size();
    if(ofs.is_open()) {
        for(auto p: result) {
            if(p.second >= 0.3) {
                ofs << seqList.list[p.first/base].header << "\t" << seqList.list[p.first%base].header 
                    << "\t" << p.second << "\n";
                pairNum ++;
            }
        }
        ofs.flush();
        ofs.close();
        std::cout << "+++ Finished LSH, " << pairNum << " possible homo pairs (including false positve) are found." << std::endl;
    } else {
        std::cerr << "!!! Error: failed to open " << outputFilename << std::endl;
    }
}