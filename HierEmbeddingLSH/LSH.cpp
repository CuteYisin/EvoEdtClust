#include "LSH.h"


DSU::DSU(int n) {
    fa = std::vector <int> (n);
    rank = std::vector <int> (n);
    for(int i = 0; i < n; i ++) {
        fa[i] = i, rank[i] = 1;
    }
}


DSU::~DSU() { }


int DSU::find(int x) {
    if(x == fa[x]) return x;
    return fa[x] = find(fa[x]);
}


void DSU::merge(int x, int y) {
    int faX = find(x), faY = find(y);
    if(faX == faY) return;
    if(rank[faX] < rank[faY]) {
        fa[faX] = faY;
    } else if(rank[faX] > rank[faY]) {
        fa[faY] = faX;
    } else {
        fa[faX] = faY;
        rank[faY] ++;
    }
}


bool DSU::isConnected(int x, int y) {
    return find(x) == find(y);
}


Fingerprint::Fingerprint(int Q = 0) {
    key = std::vector <int> (Q);
}


Fingerprint::~Fingerprint() { }


bool Fingerprint::operator<(const Fingerprint& o) const {
    return key < o.key;
}


PStableLSH::PStableLSH(const SequenceClusterNode& seqClusterNode, const GappedKmerEmbedding& gkmerEmbedding):
    seqClusterNode(seqClusterNode),
    gkmerEmbedding(gkmerEmbedding),
    dsu(seqClusterNode.n) {
        std::cout << "=== Partitioning this set use p-stable LSH " << std::endl;
        std::cout << "+++ use " << Q << " hash values as a key" << std::endl;
        sigma = getOptimalSigma();
        std::cout << "+++ Due to " << seqClusterNode.n << " sequences to cluster, optimal sigma is " << sigma << std::endl;
        T = getOptimalT();
        std::cout << "+++ Then the optimal T is " << T << std::endl;
}


PStableLSH::~PStableLSH() { }


double PStableLSH::cauchySample() {
    std::uniform_real_distribution<double> uniDist(0.0, 1.0);
    std::random_device rd;
    std::mt19937 gen(rd());
    return std::tan(std::acos(-1.0)*(uniDist(gen) - 0.5));
}


double PStableLSH::__f(double x) {
    return 2.0 * std::atan(x) / std::acos(-1.0) - (1.0 / (std::acos(-1.0)*x)) * std::log(1.0 + x * x);
}


double PStableLSH::getOptimalSigma() {
    targetFP = 0.01 / seqClusterNode.n;
    targetP2 = std::pow(targetFP, 1.0 / Q);
    double l = 1e-6, r = 100.0, mid;
    while(l - r < - 1e-4) {
        mid = (r + l) / 2.0;
        if(__f(mid) - targetP2 < - 1e-4) {
            l = mid + 1e-6;
        } else {
            r = mid - 1e-6;
        }
    }
    return mid;
}


int PStableLSH::getOptimalT() {
    return 10;
}


void PStableLSH::work() {
    for(int iter = 1; iter <= T; ++ iter) {
        std::cout << "--- Now processing iteration " << iter << std::endl;
        std::vector <Fingerprint> table = std::vector <Fingerprint> (seqClusterNode.n);
        for(int i = 0; i < seqClusterNode.n; i ++) {
            table[i] = Fingerprint(Q);
        }

        for(int k = 0; k < Q; k ++) {
            std::uniform_real_distribution<double> sigmaDist(0.0, sigma);
            std::random_device rd;
            std::mt19937 gen(rd());
            double b = sigmaDist(gen);
            std::unordered_map <HASHED_KMER, double> a;

            for(int i = 0; i < seqClusterNode.n; i ++) {
                double h = 0.0;
                auto& vec = gkmerEmbedding.profile[i];
                for(auto& v: vec) {
                    auto gkmer = v.first;
                    auto cnt = v.second;
                    if(a.count(gkmer) == 0) {
                        a[gkmer] = cauchySample();
                    }
                    h += a[gkmer] * cnt;
                }
                h /= gkmerEmbedding.nGKmer[i] + gkmerEmbedding.minNGKmer;
                table[i].key[k] = (int)((h + b) / sigma);
            }
        }

        

        for(int i = 0; i < seqClusterNode.n; i ++) {
            int seqId = seqClusterNode.seqSet[i];
            std::cout << seqId << ": " << seqClusterNode.seqList.list[seqId] << std::endl;
            for(int k = 0; k < Q; k ++) {
                std::cout << table[i].key[k] << ", ";
            }
            std::cout << std::endl;
        }
        

        std::map <Fingerprint, int> bucket;
        for(int i = 0; i < seqClusterNode.n; i ++) {
            if(bucket.count(table[i]) == 0) {
                bucket[table[i]] = dsu.find(i);
            } else {
                dsu.merge(bucket[table[i]], i);
            }
        }
    }

    finalSubset = std::unordered_map <int, std::vector<int> > ();
    for(int i = 0; i < seqClusterNode.n; i ++) {
        if(finalSubset.count(dsu.fa[i]) == 0) {
            finalSubset[dsu.fa[i]] = std::vector <int> ();
        }
        finalSubset[dsu.fa[i]].emplace_back(seqClusterNode.seqSet[i]);
    }
}