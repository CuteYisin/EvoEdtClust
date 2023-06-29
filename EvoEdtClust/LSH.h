#ifndef __LSH_H__
#define __LSH_H__


#include <list>
#include <cmath>
#include <random>
#include <ctime>

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


typedef unsigned int HASHED_KMER; 


class PStableLSH {
    private:
        double genCauchy();
        double genUniInSigma();

        inline HASHED_KMER animoAcidHash(const char&);

        double __f(double);
        double getOptimalSigma();
        int getOptimalT();


    public:
        const ClusterNode& node;
        const GappedKmerEmbedding& gke;

        double targetFP, targetP2, targetP1;
        double sigma;
        int Q, T;

        DSU dsu;
        std::unordered_map <int, std::vector<int> > subIdList;

        std::unordered_map <HASHED_KMER, std::vector <double> > a;
        std::vector <double> b;

        PStableLSH(const ClusterNode&, const GappedKmerEmbedding&);
        ~PStableLSH();

        void scanPars();
        void updatePars();

        void work();
};




#endif