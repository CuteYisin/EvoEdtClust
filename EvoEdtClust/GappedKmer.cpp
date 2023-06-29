#include "GappedKmer.h"


double ParameterGenerator::computeSimFromWK(const int& w, const int& k, const double& avgL) {
    return 1.0 - (1.0-w*(k-1)/avgL/k) / k;
}


//w as small as possible
//k as huge as possible
//void ParameterGenerator::updateWKFromSim(const double& sim, int& w, int& k, const double& avgL) {
//    k = 2;
//    w = avgL * k/(k-1) * (1.0-(1.0-sim)*k);
//    while(k <= 6 && (w < k || w > avgL/2)) {
//        k ++;
//        w = avgL * k/(k-1) * (1.0-(1.0-sim)*k);
//    }
//    if(k > 6) {
//        std::cerr << "!!! Cannot find good w and k with similarity " << sim << std::endl;
//        exit(-1);
//    }
//}

void ParameterGenerator::updateWKFromSim(const double& sim, int& w, int& k, const double& avgL) {
    for(k = 4; k >= 2; k --) {
        w = avgL * k / (k - 1) * (1.0 - (1.0 - sim) * k);
        if(w >= k && w < avgL) {
            break;
        }
    }
    if(k < 2) {
        std::cerr << "!!! Cannot find good w and k with similarity " << sim << std::endl;
        exit(-1);
    }
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
    for(int i = 0; i < node.n; ++ i) {
        const int id = node.idList[i];
        const std::string_view s = node.seqList.data[id];
        nGKmer[i] = calcNGKmer(s.length());
        maxNGKmer = std::max(maxNGKmer, nGKmer[i]);
        minNGKmer = std::min(minNGKmer, nGKmer[i]);
        avgNGKmer += nGKmer[i];
    }
    avgNGKmer /= node.n;
}