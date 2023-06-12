#ifndef __Norm_DistanceSpaceEmbedding_H__
#define __Norm_DistanceSpaceEmbedding_H__


#include <cmath>
#include "GappedKmerScan.h"


class Norm {
    public:
        const GappedKmerScan& scanner;

        Norm(const GappedKmerScan&);
        ~Norm();
};


class L1: public Norm {
    public:
        L1(const GappedKmerScan&);
        ~L1();

        void calculate(SequencePair&);
};


class L2: public Norm {
    public:
        L2(const GappedKmerScan&);
        ~L2();

        void calculate(SequencePair&);
};


class Cosine: public Norm {
    public:
        Cosine(const GappedKmerScan&);
        ~Cosine();

        void calculate(SequencePair&);
};


class Jaccard: public Norm {
    public:
        Jaccard(const GappedKmerScan&);
        ~Jaccard();

        void calculate(SequencePair&);
};

#endif
