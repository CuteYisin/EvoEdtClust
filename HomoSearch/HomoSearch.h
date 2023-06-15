#ifndef __HomoSearch_H__
#define __HomoSearch_H__


#include <cmath>
#include <random>
#include <numeric>
#include <cassert>
#include <utility>
#include <unordered_map>
#include <queue>
#include "other_unordered_map.h"
#include "GappedKmerScan.h"
#include "UnionFind.h"


typedef std::pair <int, int> PAIR;
typedef unsigned long long PAIR_TYPE;
typedef std::unordered_map <PAIR_TYPE, double> PAIR_MAP;


class HomoPairs {
    public:
        PAIR_MAP pairs;

        HomoPairs();
        ~HomoPairs();
        void addPair(const int&, const int&, const PAIR_TYPE&, double);
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
        void divideBuckets(const int&, const std::vector <int>&);

    public:
        std::vector <PAIR> result;

        LSH(const SequenceList&, const GappedKmerScan&, const int&, const double&);
        ~LSH();
        
        void work();
        void dumpToFile(const std::string&);
};


#endif
