#include "GappedKmerScan.h"


GappedKmerScan::GappedKmerScan(const int& windowSize, const int& kmerLength): 
    windowSize(windowSize),
    kmerLength(kmerLength) {
        std::cout << "\n=== Scan sequences into gapped kmer profile ===" << std::endl;
        std::cout << "+++ Initialized window size " << windowSize << ", kmer length " << kmerLength << std::endl;
        //when scanning, the current animo acid must be included, so we can only store the next k-1 positions
        otherIndices = std::vector <std::vector <int> > ();
        generateIndices();

        profile = std::vector <PROFILE_TYPE> ();
    };


GappedKmerScan::~GappedKmerScan() { }


void GappedKmerScan::generateIndices() {
    std::vector <int> indices;
    for(int i = 1; i < windowSize; ++ i) {
        dfsIndices(i, indices);
    }
}


void GappedKmerScan::dfsIndices(int now, std::vector <int>& indices) {
    indices.emplace_back(now);
    if((int)indices.size() == kmerLength - 1) {
        otherIndices.emplace_back(indices);
        indices.pop_back();
        return;
    }
    for(int j = now+1; j < windowSize; ++ j) {
        dfsIndices(j, indices);
    }
    indices.pop_back();
}


void GappedKmerScan::scan(const SequenceList& seqList) {
    int process = 0;
    for(const Sequence& s: seqList.list) {
        addSequenceProfile(s.sequence);

        process ++;
        if(process % 10000 == 0) {
            std::cout << "+++ " << process << " sequences scanned ..." << std::endl;
        }
    }
}


void GappedKmerScan::addSequenceProfile(const std::string& sequence) {
    PROFILE_TYPE temp;
    int l = sequence.length();
    for(int p = 0; p < l; ++ p) {
        for(auto& indices: otherIndices) {
            if(p + indices.back() < l) { // The last window can only generate one gapped kmer
                HASH_TYPE hashValue(1, sequence[p]);
                for(auto& ids: indices) {
                    hashValue.push_back(sequence[p+ids]);
                }
                temp[hashValue] ++;
            }
        }
        if(kmerLength == 1) temp[HASH_TYPE(1, sequence[p])] ++;
    }
    profile.emplace_back(temp);
}
