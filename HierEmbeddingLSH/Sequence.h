#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__


#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <queue>

//#include <cstring>
//#include <algorithm>
//#include <sys/mman.h>
//#include <sys/stat.h>
//#include <fcntl.h>


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
        void dumpIndexedHeader(const std::string&);
};


class LiteSequenceList {
    public:
        std::vector <std::string> list;

        LiteSequenceList();
        ~LiteSequenceList();
        void loadFromFasta_dumpIndexedHeader(const std::string&, const std::string&);
};


class SequenceClusterNode {
    public:
        const LiteSequenceList& seqList;
        std::vector <int> seqSet;
        int n;
        double avgL;

        int level;
        double expectedSimLowBound;

        SequenceClusterNode(const LiteSequenceList&, const std::vector <int>&, int, double);
        ~SequenceClusterNode();
        void show();
};


class SequenceClusterTree {
    public:
        std::queue <SequenceClusterNode> treeInBFS;
        //std::priority_queue;

        SequenceClusterTree(const LiteSequenceList&);
        ~SequenceClusterTree();
};


#endif