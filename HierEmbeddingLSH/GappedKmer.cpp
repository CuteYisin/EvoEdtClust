#include "GappedKmer.h"


double ParameterGenerator::computeSimFromWK(const int& w, const int& k, const double& avgL) {
    return 1.0 - (1.0-w*(k-1)/avgL/k) / k;
}


void ParameterGenerator::updateWKFromSim(const double& sim, int& w, int& k, const double& avgL) {
    k = 2;
    w = avgL * k/(k-1) * (1.0-(1.0-sim)*k);
    while(k <= 4 && (w < k || w > avgL/2)) {
        k ++;
        w = avgL * k/(k-1) * (1.0-(1.0-sim)*k);
    }
    if(k > 4) {
        std::cerr << "!!! Cannot find good w and k with similarity " << sim << std::endl;
    }
}


GappedKmerEmbedding::GappedKmerEmbedding(int w, int k): w(w), k(k) {
    std::cout << "=== Scan sequences into gapped kmer profile " << std::endl;

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


inline HASHED_KMER GappedKmerEmbedding::animoAcidHash(const char& aa) {
    return aa - 'A';
}


void GappedKmerEmbedding::scan(const SequenceClusterNode& seqClusterNode) {
    profile = std::vector <KMER_PROFILE> (seqClusterNode.n);
    nGKmer = std::vector <int> (seqClusterNode.n);
    maxNGKmer = 0, minNGKmer = 0x3f3f3f3f, avgNGKmer = 0.0;
    for(int i = 0; i < seqClusterNode.n; ++ i) {
        profile[i] = KMER_PROFILE();
        nGKmer[i] = 0;

        const int& seqId = seqClusterNode.seqSet[i];
        const std::string_view& s = seqClusterNode.seqList.list[seqId];
        const int& l =  s.length();
        for(int p = 0; p < l; p ++) {
            for(const auto& gkmerMask: otherIndices) {
                if(p + gkmerMask.back() < l ) {
                    HASHED_KMER hashValue = animoAcidHash(s[p]);
                    for(const auto& j: gkmerMask) {
                        hashValue <<= 5;
                        hashValue += animoAcidHash(s[p+j]);
                    }
                    profile[i][hashValue] ++;
                    nGKmer[i] ++;
                }
            }
        }
        maxNGKmer = std::max(maxNGKmer, nGKmer[i]);
        minNGKmer = std::min(minNGKmer, nGKmer[i]);
        avgNGKmer += nGKmer[i];
    }
    avgNGKmer /= seqClusterNode.n;

    /*
    for(int i = 0; i < seqClusterNode.n; ++ i) {
        const int& seqId = seqClusterNode.seqSet[i];
        const std::string& s = seqClusterNode.seqList.list[seqId];
        std::cout << s << std::endl;
        int tt = 0;
        for(auto& x: profile[i]) {
            std::cout << x.first << "," << x.second << "; \t";
            tt += x.second;
        }
        std::cout << std::endl;
        std::cout << "--- " << tt << " vs " << nGKmer[i] << " vs " << calcN(s.length()) << std::endl;
    }
    */
}


int GappedKmerEmbedding::calcN(int l) {
    int s = 1;
    for(int r = 1; r <= k-1 ; r ++) {
        s *= w-r;
        s /= r;
    }
    return (l-w)*s + s * w / k;
}