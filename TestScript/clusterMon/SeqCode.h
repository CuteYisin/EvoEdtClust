#ifndef SEQCODE_CLUSTERMON_H_
#define SEQCODE_CLUSTERMON_H_

#include <string>
#include <unordered_map>
#include <fstream>
#include <cstring>
#include <iostream>

typedef unsigned int ID_TYPE;

class SeqCode {
    public:
        std::unordered_map <std::string, ID_TYPE> seqID;

        SeqCode();
        ~SeqCode();

        void code(const std::string&);
};

#endif
