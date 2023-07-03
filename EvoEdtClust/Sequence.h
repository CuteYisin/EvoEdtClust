#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__


#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <queue>


class SequenceList {
    public:
        std::vector <std::string> data;

        SequenceList();
        ~SequenceList();
        void loadFromFasta_dumpIndexedHeader(const std::string&, const std::string&);
};


class ClusterNode {
    public:
        const SequenceList& seqList;
        std::vector <int> idList;
        int level;
        double expectedSimLowBound;

        int n;
        double avgL;

        ClusterNode(const SequenceList&, const std::vector <int>&, int, double);
        ~ClusterNode();

        //bool operator< (const ClusterNode&) const; // return n > other.n;

        void show();
};


class ClusterTree {
    public:
        std::queue <ClusterNode> bfsOrder;
        //std::priority_queue <ClusterNode> bfsOrder;

        ClusterTree(const SequenceList&);
        ~ClusterTree();
};


#endif