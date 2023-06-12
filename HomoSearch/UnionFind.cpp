#include "UnionFind.h"

UnionFind::UnionFind(const int& pointNumber) {
    parent = std::vector <int> (pointNumber);
    for(int i = 0; i < pointNumber; i++) {
        parent[i] = i;
    }

    rank = std::vector <int> (pointNumber, 0);
    count = pointNumber;
}


UnionFind::~UnionFind() {}


int UnionFind::find(const int& x) {
    if(parent[x] != x) {
        parent[x] = find(parent[x]);
    }
    return parent[x];
}


bool UnionFind::connected(const int& x, const int& y) {
    return find(x) == find(y);
}


void UnionFind::merge(const int& x, const int& y) {
    int parentX = find(x);
    int parentY = find(y);
    if(parentX == parentY) return;
    if(rank[parentX] < rank[parentY]) {
        parent[parentX] = parentY;
    } else if (rank[parentX] > rank[parentY]) {
        parent[parentY] = parentX;
    } else {
        parent[parentY] = parentX;
        rank[parentX] ++;
    }
    count --;
}