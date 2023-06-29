#ifndef __GappedKmer_H__
#define __GappedKmer_H__


#include <unordered_map>
#include "Sequence.h"


typedef unsigned int HASHED_KMER;
typedef std::unordered_map <HASHED_KMER, int> KMER_PROFILE;


class ParameterGenerator {
    public:
        static double computeSimFromWK(const int&, const int&, const double&);
        static void updateWKFromSim(const double&, int&, int&, const double&);
};


class GappedKmerEmbedding {
    private:
        std::vector <std::vector <int> > otherIndices;

        void dfsIndices(int, std::vector <int>&);
        inline HASHED_KMER animoAcidHash(const char&);

    public:
        const int w, k;
        std::vector <KMER_PROFILE> profile;
        std::vector <int> nGKmer;
        int maxNGKmer, minNGKmer;
        double avgNGKmer;

        GappedKmerEmbedding(int, int);
        ~GappedKmerEmbedding();

        void scan(const SequenceClusterNode&);
        int calcN(int);
};


#endif