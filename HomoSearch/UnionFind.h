#ifndef __UnionFind_H__
#define __UnionFind_H__

#include <vector>

class UnionFind {
    public:
        std::vector<int> parent;
        std::vector<int> rank;
        int count;

        UnionFind(const int&);
        ~UnionFind();

        int find(const int&);
        bool connected(const int&, const int&);
        void merge(const int&, const int&);
};

#endif