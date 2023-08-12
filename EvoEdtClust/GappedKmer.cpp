#include <cmath>
#include "GappedKmer.h"

//w-k (w,k)
std::vector <int> ParameterGenerator::w = {2, 5};
std::vector <int> ParameterGenerator::k = {2, 3};


double ParameterGenerator::computeSimFromWK(const int& w, const int& k, const double& avgL) {
    return 1.0 - (1.0-w*(k-1)/avgL/k) / k;
}


GappedKmerEmbedding::GappedKmerEmbedding(const ClusterNode& node, int w, int k): 
    node(node), w(w), k(k) {
    std::cout << "=== Construct gapped kmer scanner " << std::endl;

    otherIndices = std::vector <std::vector <int> > ();
    std::vector <int> indices;
    for(int i = 1; i < w; i ++) {
        dfsIndices(i, indices);
    }
    std::cout << "+++ " << otherIndices.size() << " gapped " << k << "-mers will be scanned in a " << w << "-length window" << std::endl;
}


GappedKmerEmbedding::~GappedKmerEmbedding() { }


void GappedKmerEmbedding::dfsIndices(int now, std::vector <int>& indices) {
    indices.emplace_back(now);
    if((int)indices.size() == k-1) {
        otherIndices.emplace_back(indices);
        indices.pop_back();
        return;
    }
    for(int j = now+1; j < w; ++ j) {
        dfsIndices(j, indices);
    }
    indices.pop_back();
}


int GappedKmerEmbedding::calcNGKmer(int l) {
    int t = 1;
    for(int r = 1; r <= k-1 ; r ++) {
        t *= w-r;
        t /= r;
    }
    return (l-w)*t + t * w / k;
}


void GappedKmerEmbedding::scan() {
    nGKmer = std::vector <int> (node.n);
    maxNGKmer = 0, minNGKmer = 0x3f3f3f3f, avgNGKmer = 0.0;
    theoreticalTotalNGKMer = std::pow(20, k);
    for(int i = 0; i < node.n; ++ i) {
        const int id = node.idList[i];
        const std::string_view s = node.seqList.data[id];
        nGKmer[i] = calcNGKmer(s.length());
        maxNGKmer = std::max(maxNGKmer, nGKmer[i]);
        minNGKmer = std::min(minNGKmer, nGKmer[i]);
        avgNGKmer += nGKmer[i];
    }
    theoreticalTotalNGKMer = std::min(theoreticalTotalNGKMer, avgNGKmer);
    avgNGKmer /= node.n;
}