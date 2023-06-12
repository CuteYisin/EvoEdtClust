#include "LSH.h"


LSH::LSH(const int& repetition, const GappedKmerScan& scanner): repetition(repetition), scanner(scanner) {
    r_opt = 0.1;
    double pi = std::acos(-1.0);
    double p2 = (2 * std::atan(4 * r_opt) / pi) - 1.0 / (pi * 4 * r_opt) * std::log(1 + 16 * r_opt * r_opt);
    k_opt = std::round(std::log(scanner.profile.size()) / (-std::log(p2)));
}


LSH::~LSH() {}


void LSH::print(const std::string& outputFilename, int mode) {
    std::ofstream ofs(outputFilename);
    int seqNumber = scanner.profile.size();
    if(ofs.is_open()) {
        if(mode == 1) {
            for(int i = 0; i < seqNumber; i++) {
                for(int j = i + 1; j < seqNumber; j++) {
                    int value = 0;
                    for(int l = 0; l < repetition; l++) {
                        int conflictNumber = 0;
                        for(int k = 0; k < k_opt; k++) {
                            if(abs(estimationValueDouble[l][i][k] - estimationValueDouble[l][j][k]) <= r_opt) {
                                conflictNumber ++;
                            }
                        }
                        if(conflictNumber >= k_opt * 0.8) {
                            value ++;
                        }
                    }
                    ofs << repetition - value << "\t";
                }
                ofs << std::endl;
            }
        } else if(mode == 2) {
            for(int i = 0; i < seqNumber; i++) {
                for(int j = i + 1; j < seqNumber; j++) {
                    int value = 0;
                    for(int k = 0; k < repetition; k++) {
                        if(estimationValueInt[k][i] != estimationValueInt[k][j]) {
                            value ++;
                        }
                    }
                    ofs << value << "\t";
                }
                ofs << std::endl;
            }
        } else {
            for(int i = 0; i < seqNumber; i++) {
                for(int j = i + 1; j < seqNumber; j++) {
                    int value = 0;
                    std::unordered_map <int, int> findMap;
                    for(int k = 0; k < repetition; k++) {
                        findMap[estimationValueInt[k][i]] ++;
                    }
                    for(int k = 0; k < repetition; k++) {
                        if(findMap[estimationValueInt[k][j]]) {
                            value ++;
                            findMap[estimationValueInt[k][j]] --;
                        }
                    }
                    ofs << repetition - value << "\t";
                }
                ofs << std::endl;
            }
        }
        ofs.flush();
        ofs.close();
    } else {
        std::cerr << "!!! Error: failed to open " << outputFilename << std::endl;
    }
}


pStable::pStable(const int& repetition, const GappedKmerScan& scanner): LSH(repetition, scanner) {
    L1Estimation = std::vector <std::vector <std::vector <double>>>
        (repetition, std::vector <std::vector <double>>(scanner.profile.size(), std::vector <double>()));
    L2Estimation = std::vector <std::vector <std::vector <double>>>
        (repetition, std::vector <std::vector <double>>(scanner.profile.size(), std::vector <double>()));
}


pStable::~pStable() {}


void pStable::work(int p) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::cauchy_distribution<double> cauchy_dist(0.0, 1.0);
    std::normal_distribution<double> normal_dist(0.0, 1.0);

    std::vector <long long> kmerTotalNumber;
    for(auto x: scanner.profile) {
        int number = 0;
        for(auto kmer: x) {
            number += kmer.second;
        }
        kmerTotalNumber.emplace_back(number);
    }

    for(int i = 0; i < repetition; i++) {
        std::vector <std::vector <double>> dot_products(scanner.profile.size(), std::vector <double>());
        for(int k = 0; k < k_opt; k ++) {
            std::unordered_map <HASH_TYPE, double> Table;
            int index = 0;
            for(auto x: scanner.profile) {
                double dot_product = 0;
                for(auto kmer: x)   {
                    if(Table.count(kmer.first) == 0) {
                        double value;
                        if(p == 1) {
                            value = cauchy_dist(gen);
                        } else {
                            value = normal_dist(gen);
                        }
                        Table[kmer.first] = value;
                    }
                    dot_product += kmer.second * Table[kmer.first];
                }
                dot_products[index].emplace_back(dot_product/(2 * kmerTotalNumber[index]));
                index ++;
            }
        }
        int sequenceNumber = scanner.profile.size();
        if(p == 1) {
            for(int x = 0; x < sequenceNumber; x ++) {
                L1Estimation[i][x] = dot_products[x];
            }
        } else {
            for(int x = 0; x < sequenceNumber; x ++) {
                L2Estimation[i][x] = dot_products[x];
            }
        }
    }
    if(p == 1) {
        estimationValueDouble = L1Estimation;
    } else {
        estimationValueDouble = L2Estimation;
    }
}


RandomProjection::RandomProjection(const int& repetition, const GappedKmerScan& scanner): LSH(repetition, scanner) {
    CosineEstimation = std::vector <std::vector <int>> (repetition, std::vector <int>());
}


RandomProjection::~RandomProjection() {}


void RandomProjection::work() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(1, 6);

    for(int i = 0; i < repetition; i++) {
        std::unordered_map <HASH_TYPE, double> Table;
        
        for(auto x: scanner.profile) {
            double dot_product = 0;
            for(auto kmer: x)   {
                if(Table.count(kmer.first) == 0) {
                    double value;
                    int randomInt = distribution(gen);
                    if(randomInt == 1) {
                        value = sqrt(3);
                    } else if (randomInt < 6) {
                        value = 0;
                    } else {
                        value = -sqrt(3);
                    }
                    Table[kmer.first] = value;
                }
                dot_product += kmer.second * Table[kmer.first];
            }
            if(dot_product >= 0) {
                CosineEstimation[i].emplace_back(1);
            } else {
                CosineEstimation[i].emplace_back(-1);
            }
        }
    }
    estimationValueInt = CosineEstimation;
}


MinHash::MinHash(const int& repetition, const GappedKmerScan& scanner): LSH(repetition, scanner) {
    JaccardEstimation = std::vector <std::vector <int>> (repetition, std::vector <int>());
}


MinHash::~MinHash() {}


int MinHash::hash_djb2(const HASH_TYPE& str) {
    int hash = 5381;
    int c, index = 0;

    while ((c = str[index++])) {
        hash = ((hash << 5) + hash) + c;
    }

    return hash;
}


void MinHash::work() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution;
    int MOD = distribution(gen);

    for(auto x: scanner.profile) {
        std::priority_queue <int, std::vector<int>, std::greater<int>> minValue;
        for(auto kmer: x) {
            HASH_TYPE str = kmer.first;
            int value = hash_djb2(str) % MOD;
            for(int i = 0; i < kmer.second; i++) {
                minValue.emplace(value);
            }
        }
        for(int i = 0; i < repetition; i++) {
            JaccardEstimation[i].emplace_back(minValue.top());
            minValue.pop();
        }
    }
    estimationValueInt = JaccardEstimation;
}