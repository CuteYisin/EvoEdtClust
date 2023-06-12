#ifndef __Sequence_H__
#define __Sequence_H__


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


#endif
