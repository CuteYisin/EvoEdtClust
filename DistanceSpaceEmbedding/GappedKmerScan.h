#ifndef __GappedKmerScan_DistanceSpaceEmbedding_H__
#define __GappedKmerScan_DistanceSpaceEmbedding_H__


#include <unordered_map>
#include "Sequence.h"


typedef std::string HASH_TYPE;
typedef std::unordered_map <HASH_TYPE, int> PROFILE_TYPE;


class GappedKmerScan {
    private:
        std::vector <std::vector <int> > otherIndices;

        void generateIndices();
        void dfsIndices(int, std::vector <int>&);

        void addSequenceProfile(const std::string&);

    public:
        const int windowSize, kmerLength;

        std::vector <PROFILE_TYPE> profile;

        GappedKmerScan(const int&, const int&);
        ~GappedKmerScan();
        void scan(const SequenceList&);
};

#endif
