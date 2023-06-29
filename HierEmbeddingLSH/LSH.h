#ifndef __LSH_H__
#define __LSH_H__


#include <list>
#include <cmath>
#include <random>
#include <map>

#include "Sequence.h"
#include "GappedKmer.h"


class DSU {
    public:
        std::vector <int> fa, rank;

        DSU(int);
        ~DSU();

        int find(int);
        void merge(int, int);
        bool isConnected(int, int);
};


class Fingerprint {
    public:
        std::vector <int> key;

        Fingerprint(int);
        ~Fingerprint();

        bool operator<(const Fingerprint&) const;
};


class PStableLSH {
    private:
        double cauchySample();
        double __f(double);
        double getOptimalSigma();
        int getOptimalT();

    public:
        const SequenceClusterNode& seqClusterNode;
        const GappedKmerEmbedding& gkmerEmbedding;

        
        int Q = 3;
        double targetFP, targetP2, targetP1;
        double sigma;
        int T;

        DSU dsu;
        std::unordered_map <int, std::vector<int> > finalSubset;

        PStableLSH(const SequenceClusterNode&, const GappedKmerEmbedding&);
        ~PStableLSH();

        void work();
};




#endif