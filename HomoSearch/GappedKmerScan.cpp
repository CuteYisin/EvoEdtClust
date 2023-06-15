#include "GappedKmerScan.h"


GappedKmerScan::GappedKmerScan(std::vector <int>& windowSize, std::vector <int>& kmerLength):
    windowSize(windowSize),
    kmerLength(kmerLength),
    levelNumber(windowSize.size()) {
        std::cout << "\n=== Scan sequences into gapped kmer profile ===" << std::endl;

        otherIndices = std::vector <std::vector <std::vector <int> > > (levelNumber, std::vector <std::vector <int> >());
        for(int i = 0; i < levelNumber; i++) {
            std::cout << "+++ Initialized window size " << windowSize[i] << ", kmer length " << kmerLength[i] << std::endl;
            //when scanning, the current animo acid must be included, so we can only store the next k-1 positions
            generateIndices(i);
        }
        profile = std::vector <std::vector <PROFILE_TYPE> > ();
    }


GappedKmerScan::~GappedKmerScan() { }


void GappedKmerScan::generateIndices(int levelNow) {
    std::vector <int> indices;
    for(int i = 1; i < windowSize[levelNow]; ++ i) {
        //Start with i and search for index combinations that meet the w constraints
        dfsIndices(levelNow, i, indices);
    }
}


void GappedKmerScan::dfsIndices(int levelNow, int now, std::vector <int>& indices) {
    indices.emplace_back(now);
    if((int)indices.size() == kmerLength[levelNow] - 1) {
        otherIndices[levelNow].emplace_back(indices);    //Store eligible index combinations
        indices.pop_back();
        return;
    }
    for(int j = now+1; j < windowSize[levelNow]; ++ j) {
        dfsIndices(levelNow, j, indices);
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
    int l = sequence.length();
    std::vector <PROFILE_TYPE> temp(levelNumber, PROFILE_TYPE());
    for(int i = 0; i < levelNumber; i++) {
        int k = kmerLength[i];
        for(int p = 0; p < l; ++ p) {
            for(auto& indices: otherIndices[i]) {
                // The last window can only generate one gapped kmer
                if(p + indices.back() < l) { 
                    HASH_TYPE hashValue(1, sequence[p]);
                    for(auto& ids: indices) {
                        hashValue.push_back(sequence[p+ids]);
                    }
                    temp[i][hashValue] ++;
                } 
            }
            if(k == 1) temp[i][HASH_TYPE(1, sequence[p])] ++;
        }
    }
    profile.emplace_back(temp);
}
