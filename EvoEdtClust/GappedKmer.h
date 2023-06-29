#ifndef __GappedKmer_H__
#define __GappedKmer_H__


#include <unordered_map>
#include "Sequence.h"


class ParameterGenerator {
    public:
        static double computeSimFromWK(const int&, const int&, const double&);
        static void updateWKFromSim(const double&, int&, int&, const double&);
};


class GappedKmerEmbedding {
    private:
        void dfsIndices(int, std::vector <int>&);

    public:
        const ClusterNode& node;
        const int w, k;

        std::vector <std::vector <int> > otherIndices;

        std::vector <int> nGKmer;
        int maxNGKmer, minNGKmer;
        double avgNGKmer;

        GappedKmerEmbedding(const ClusterNode&, int, int);
        ~GappedKmerEmbedding();

        int calcNGKmer(int);

        void scan();
};


#endif