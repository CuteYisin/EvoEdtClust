#ifndef __Sequence_DistanceSpaceEmbedding_H__
#define __Sequence_DistanceSpaceEmbedding_H__


#include <fstream>
#include <iostream>
#include <string>
#include <vector>


class Sequence {
    public:
        std::string header, sequence;
        int index;

        Sequence(const std::string&, const std::string&, const int&);
        ~Sequence();
};


class SequenceList {
    public:
        std::vector <Sequence> list;

        SequenceList();
        ~SequenceList();
        void loadFromFasta(const std::string&);
};


class SequencePair {
    public:
        const int pairNumber;
        std::vector <std::vector <double> > editDistances, L1Distances, L2Distances, cosineDistances, jaccardDistances;

        SequencePair(const SequenceList&);
        ~SequencePair();
        void addPair(int&, int&, double&, std::vector <std::vector <double> >&);
        void dumpToFile(const std::string&, std::vector <std::vector <double> >&);
};


#endif