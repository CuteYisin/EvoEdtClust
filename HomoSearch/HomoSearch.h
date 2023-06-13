#ifndef __HomoSearch_H__
#define __HomoSearch_H__


#include <cmath>
#include <random>
#include <numeric>
#include <cassert>
#include <utility>
#include <unordered_map>
#include <queue>
#include "GappedKmerScan.h"
#include "UnionFind.h"


typedef unsigned long long PAIR_TYPE;
typedef std::unordered_map <PAIR_TYPE, double> PAIR_MAP;


class HomoPairs {
    public:
        PAIR_MAP pairs;
        PAIR_TYPE base;

        HomoPairs(PAIR_TYPE&);
        ~HomoPairs();
        void addPair(const int&, const int&, double);
};


class HomoSearch {
    protected:
        const SequenceList& seqList;

    public:
        HomoPairs homoPairs();

        HomoSearch(const SequenceList&);
        ~HomoSearch();
};


class LSH: public HomoSearch {
    private:
        const GappedKmerScan& scanner;
        const int mode;
        const int pstable_hashTable_k, pstable_hashTable_L, minHash_hashTable_k, minHash_hashRable_L;
        const double similarityThreshold, projectionWidth;

        const int table_size;
        const double location, scale, table_step;

        std::vector<double> cauchyTable;
        std::vector< std::vector <int> > combinations;

        void makeCauchyTable();
        void generateCombinations();
        int hashFunction(const HASH_TYPE&);

        int getFingerprint_pStable(const int&, const std::vector <PROFILE_TYPE>&, 
            std::unordered_map <HASH_TYPE, double>&, const int&, const double&);
        std::vector <int> getFingerprint_minHash(const int&, const std::vector <PROFILE_TYPE>&, const int&, const int&);
        PAIR_MAP divideBuckets(const int&, const std::vector <int>&, const PAIR_MAP&);

    public:
        PAIR_MAP result;

        LSH(const SequenceList&, const GappedKmerScan&, const int&, const double&);
        ~LSH();
        
        void work();
        void dumpToFile(const std::string&);
};


#endif
