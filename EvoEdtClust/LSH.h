#ifndef __LSH_H__
#define __LSH_H__


#include <list>
#include <cmath>
#include <random>
#include <ctime>
#include <map>
#include <set>
#include <algorithm>

#include "Sequence.h"
#include "GappedKmer.h"


class DSU {
    public:
        std::vector <int> fa, rank;

        DSU(int);
        ~DSU();

        int find(int);
        void merge(int, int);
};


class Fingerprint {
    public:
        std::vector <int> key;

        Fingerprint(int);
        ~Fingerprint();

        bool operator==(const Fingerprint&) const;
};


typedef unsigned long long HASHED_KMER; 


class PStableLSH {
    private:
        static int rdSeed;
        int nextRdSeed();
        double genCauchy();
        double genUniInSigma();

        inline HASHED_KMER animoAcidHash(const char&);

        double __f(double);
        double __rho(double, double, double);
        double getQ5(int, int, double);
        double getQ80(int, int, double);
        double getOptimalSigma();
        int getOptimalQ(double);
        int getOptimalT(double, int);


    public:
        const ClusterNode& node;
        const GappedKmerEmbedding& gke;

        double similarityPairEstimation;
        double targetP1, targetP2, Q5, Q80, targetRho;
        double sigma, gamma;
        int Q, T;
        
        DSU dsu;
        std::unordered_map <int, std::vector<int> > subIdList;
        // std::map <int, std::vector<int> > subIdList;

        std::vector <int> maxT;
        std::unordered_map <HASHED_KMER, std::vector <double> > a;
        std::vector <double> b;

        PStableLSH(const ClusterNode&, const GappedKmerEmbedding&);
        ~PStableLSH();

        double editDistanceCalculate(const std::string&, const std::string&);
        double sample();

        void scanPars();
        void updatePars();

        void work();
};




#endif