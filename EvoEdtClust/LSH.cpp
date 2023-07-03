#include <omp.h>
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


Fingerprint::Fingerprint(int Q = 0) {
    key = std::vector <int> (Q);
}


Fingerprint::~Fingerprint() { }


bool Fingerprint::operator==(const Fingerprint& o) const {
    return key == o.key;
}


// Hash for a vector, defined for unordered_map <Fingerprint, int> bucket;
struct FingerprintHasher {
    std::size_t operator() (const Fingerprint& x) const {
        std::size_t hash = x.key.size();
        for(auto i: x.key) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};


PStableLSH::PStableLSH(const ClusterNode& node, const GappedKmerEmbedding& gke):
    node(node), gke(gke), Q(3), T(5), dsu(node.n) {
        std::cout << "=== Partitioning this set use p-stable LSH " << std::endl;
        std::cout << "+++ use " << Q << " hash values as a key, repeat " << T << " times" << std::endl;
        sigma = getOptimalSigma();
        std::cout << "+++ Due to " << node.n << " sequences to cluster, optimal sigma is " << sigma << std::endl;
}


PStableLSH::~PStableLSH() { }


int PStableLSH::rdSeed = 1;


int PStableLSH::nextRdSeed() {
    rdSeed = rdSeed * 1664525 + 1013904223;
    return rdSeed;
}


double PStableLSH::__f(double x) {
    return 2.0 / std::acos(-1.0) * std::atan(x) - (1.0 / (std::acos(-1.0) * x)) * std::log(1.0 + x * x);
}


// hard part ***
double PStableLSH::getOptimalSigma() {
    //targetFP = 0.1 / node.n;
    //targetFP = std::pow(1.0 / node.n, 1.1);
    targetFP = 1.0 / (double)(node.n-1) / std::log(node.avgL); // if FP related with AvgL; Prior?
    //targetFP = 0.05 / (double)(node.n-1);
    //targetP2 = std::pow(targetFP, 1.0 / Q);
    targetP2 = std::pow(1.0-std::pow(1.0-targetFP, 1.0/T), 1.0/Q);
    //std::cerr << "~~~ oldp2 = " << oldP2 << ", targetp2 = " << targetP2 << std::endl;
    double l = 1e-6, r = 100.0, mid;
    while(l - r < - 1e-4) {
        mid = (r + l) / 2.0;
        if(__f(mid) - targetP2 < - 1e-4) {
            l = mid + 1e-6;
        } else {
            r = mid - 1e-6;
        }
    }
    //return mid * 2.0 * gke.calcNGKmer(median_l); 
    return mid * 2.0 * gke.avgNGKmer;
}


double PStableLSH::genCauchy() {
    std::uniform_real_distribution<double> uniDist(0.0, 1.0);
    //std::random_device rd;
    //std::mt19937 gen(rd());
    std::mt19937 gen(nextRdSeed());
    return std::tan(std::acos(-1.0)*(uniDist(gen) - 0.5));
}


double PStableLSH::genUniInSigma() {
    std::uniform_real_distribution<double> uniDist(0.0, sigma);
    //std::random_device rd;
    //std::mt19937 gen(rd());
    std::mt19937 gen(nextRdSeed());
    return uniDist(gen);
}


inline HASHED_KMER PStableLSH::animoAcidHash(const char& aa) {
    return aa - 'A';
}


void PStableLSH::scanPars() {
    a = std::unordered_map <HASHED_KMER, std::vector<double> > ();

    for(int i = 0; i < node.n; ++ i) {
        auto& id = node.idList[i];
        auto& s = node.seqList.data[id];
        const int& l = s.length();

        for(int p = 0; p < l; ++ p) {
            for(auto& mask: gke.otherIndices) {
                if(p + mask.back() < l) {
                    HASHED_KMER gkmer = animoAcidHash(s[p]);
                    for(auto& j: mask) {
                        gkmer <<= 5;
                        gkmer ^= animoAcidHash(s[p+j]);
                    }
                    if(a.count(gkmer) == 0) {
                        a[gkmer] = std::vector <double> (Q);
                        for(int k = 0; k < Q; ++ k) {
                            a[gkmer][k] = genCauchy();
                        }
                    }
                }
            }
        }
    }

    b = std::vector <double> (Q);
    for(int k = 0; k < Q; ++ k) {
        b[k] = genUniInSigma();
        //std::cerr << "~~~ test seed b[" << k << "]=" << b[k] << std::endl;
    }
}


void PStableLSH::updatePars() {
    for(auto& a_gkmer: a) {
        for(int k = 0; k < Q; ++ k) {
            a_gkmer.second[k] = genCauchy();
        }
    }
    for(int k = 0; k < Q; ++ k) {
        b[k] = genUniInSigma();
        //std::cerr << "~~~ test seed b[" << k << "]=" << b[k] << std::endl;
    }
}


void PStableLSH::work() {

    for(int iter = 1; iter <= T; ++ iter) {
        //std::cout << "--- Now processing iteration " << iter << std::endl;

        std::vector <Fingerprint> table = std::vector <Fingerprint> (node.n);
        for(int i = 0; i < node.n; ++ i) {
            table[i] = Fingerprint(Q);
        }

        updatePars();

        #pragma omp parallel for num_threads(16)
        for(int i = 0; i < node.n; ++ i) {
            std::vector <double> h(Q);
            for(int k = 0; k < Q; ++ k) {
                h[k] = 0.0;
            }

            auto& id = node.idList[i];
            auto& s = node.seqList.data[id];
            const int& l = s.length();

            for(int p = 0; p < l; ++ p) {
                for(auto& mask: gke.otherIndices) {
                    if(p + mask.back() < l) {
                        HASHED_KMER gkmer = animoAcidHash(s[p]);
                        for (auto& j: mask) {
                            gkmer <<= 5;
                            gkmer ^= animoAcidHash(s[p+j]);
                        }
                        for(int k = 0; k < Q; ++ k) {
                            h[k] += a[gkmer][k];
                        }
                    }
                }
            }

            for(int k = 0; k < Q; ++ k) {
                //h[k] /= gke.nGKmer[i] + gke.minNGKmer;
                //h[k] /= 2 * gke.nGKmer[i];
                table[i].key[k] = (int)((h[k] + b[k]) / sigma);
            }
        }


       // for(int i = 0; i < node.n; ++ i) {
       //     int id = node.idList[i];
       //     std::cout << id << ": " << node.seqList.data[id] << std::endl;
       //     for(int k = 0; k < Q; k ++) {
       //        std::cout << table[i].key[k] << ", ";
       //     }
       //     std::cout << std::endl;
       // }
        
        //Merge dsu if they shared the same fingerprint
        std::unordered_map <Fingerprint, int, FingerprintHasher> bucket;
        for(int i = 0; i < node.n; i ++) {
            if(bucket.count(table[i]) == 0) {
                bucket[table[i]] = dsu.find(i);
            } else {
                dsu.merge(bucket[table[i]], i);
            }
        }

        std::unordered_map <int, std::vector<int> > glanceDSU;
        for(int i = 0; i < node.n; ++ i) {
            if(glanceDSU.count(dsu.find(i)) == 0) {
                glanceDSU[dsu.find(i)] = std::vector <int> ();
            }
            glanceDSU[dsu.find(i)].emplace_back(node.idList[i]);
        }


       // std::cout << "~~~ Glance DSU" << std::endl;
       // for(auto v: glanceDSU) {
       //     std::cout << v.first << ": ";
       //     for(auto x: v.second) {
       //         std::cout << x << ", ";
       //     }
       //     std::cout << std::endl;
       // }
        
    }

    //Scan the dsu, and parse it to diffent sequence idlist
    subIdList = std::unordered_map <int, std::vector<int> > ();
    for(int i = 0; i < node.n; ++ i) {
        if(subIdList.count(dsu.find(i)) == 0) {
            subIdList[dsu.find(i)] = std::vector <int> ();
        }
        subIdList[dsu.find(i)].emplace_back(node.idList[i]);
    }
}

