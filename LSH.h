#ifndef __LSH_DistanceSpaceEmbedding_H__
#define __LSH_DistanceSpaceEmbedding_H__

#include <random>
#include <queue>
#include <unordered_map>
#include <cmath>
#include "GappedKmerScan.h"

class LSH{
    public:
        const int& repetition;
        const GappedKmerScan& scanner;
        double r_opt;
        int k_opt;
        std::vector <std::vector <std::vector <double>>> estimationValueDouble;
        std::vector <std::vector <int>> estimationValueInt;

        LSH(const int&, const GappedKmerScan&);
        ~LSH();

        void print(const std::string&, int);
};


class pStable: public LSH {
    public:
        std::vector <std::vector <std::vector <double>>> L1Estimation, L2Estimation;

        pStable(const int&, const GappedKmerScan&);
        ~pStable();

        void work(int p);
};


class RandomProjection: public LSH {
    public:
        std::vector <std::vector <int>> CosineEstimation;

        RandomProjection(const int&, const GappedKmerScan&);
        ~RandomProjection();

        void work();

};


class MinHash: public LSH {
    public:
        std::vector <std::vector <int>> JaccardEstimation;

        MinHash(const int&, const GappedKmerScan&);
        ~MinHash();

        int hash_djb2(const HASH_TYPE&);
        void work();
};

#endif