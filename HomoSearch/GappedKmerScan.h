#ifndef __GappedKmerScan_H__
#define __GappedKmerScan_H__


#include <unordered_map>
#include "Sequence.h"


typedef std::string HASH_TYPE;
typedef std::unordered_map <HASH_TYPE, int> PROFILE_TYPE;


class GappedKmerScan {
    private:
        std::vector <std::vector <std::vector <int> > > otherIndices;

        void generateIndices(int);
        void dfsIndices(int, int, std::vector <int>&);

        void addSequenceProfile(const std::string&);

    public:
        const std::vector <int> windowSize, kmerLength;
        const int levelNumber;

        //std::vector <PROFILE_TYPE> profile;
        std::vector< std::vector <PROFILE_TYPE> > profile;

        GappedKmerScan(std::vector <int>&, std::vector <int>&);
        ~GappedKmerScan();
        void scan(const SequenceList&);
};

#endif
