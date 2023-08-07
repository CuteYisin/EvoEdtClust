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


PStableLSH::PStableLSH(const ClusterNode& node, const GappedKmerEmbedding& gke, const double& sim):
    node(node), gke(gke), sim(sim), gamma(1.0), dsu(node.n), 
    rho({0.75, 0.9, 0.9, 0.9}),
    targetFN({0.01, 0.01, 0.05, 0.1}) {
        std::cout << "=== Partitioning this set use p-stable LSH " << std::endl;
        Q1 = gke.avgNGKmer * 2.0 * getQ1(gke.w, gke.k, node.avgL);
        Q5 = gke.avgNGKmer * 2.0 * getQ5(gke.w, gke.k, node.avgL);
        Q20 = gke.avgNGKmer * 2.0 * getQ20(gke.w, gke.k, node.avgL);
        Q80 = gke.avgNGKmer * 2.0 * getQ80(gke.w, gke.k, node.avgL);
        Q95 = gke.avgNGKmer * 2.0 * getQ95(gke.w, gke.k, node.avgL);
        std::cout << "--- Q1 is " << Q1 << ", Q5 is " << Q5 << ", Q20 is " << Q20 
            << ", Q80 is " << Q80 << ", and Q95 is " << Q95 << std::endl;
        sigma = getOptimalSigma();
        Q = getOptimalQ(sigma);
        T = getOptimalT(sigma, Q);
        std::cout << "+++ use " << Q << " hash values as a key, repeat " << T << " times" << std::endl;
        std::cout << "+++ Due to " << node.n << " sequences to cluster, optimal sigma is " << sigma << std::endl;
}


PStableLSH::~PStableLSH() { }


int PStableLSH::rdSeed = 1;


int PStableLSH::nextRdSeed() {
    rdSeed = rdSeed * 1664525 + 1013904223;
    return rdSeed;
}


double PStableLSH::__f(double x) {
    return 2.0 / std::acos(-1.0) * std::atan(x) - 
        (1.0 / (std::acos(-1.0) * x)) * std::log(1.0 + (x * x));
}


double PStableLSH::__rho(double sigma, double r1, double r2) {
    return std::log(1.0 / __f(sigma / r1 / gamma)) / std::log(1.0 / __f(sigma / r2 / gamma));
}


// hard part ***
double PStableLSH::getOptimalSigma() {
    double targetRho =rho[gke.k - 1];
    // similarButUnequal = 1 / (node.n * node.n * 0.01);
    // targetP1 = std::pow(1.0-std::pow(similarButUnequal, 1.0/T), 1.0/Q);

    // UnsimilarButEqual = 1 / (node.n * node.n * 0.83);
    // targetP2 = std::pow(1.0-std::pow(1.0-UnsimilarButEqual, 1.0/T), 1.0/Q);
    //targetFP = node.avgL * node.avgL / (std::pow(node.n, 2) * sim);
    //targetFP = (1 - sim) / (double)(node.n-1) / std::log(node.avgL); // if FP related with AvgL; Prior?
    //targetFP = 0.05 / (double)(node.n-1);
    //targetFP = 1.0 / (double)(node.n-1) / (1.0 - std::pow(1.0 / 20.0, node.avgL));
    //targetP2 = std::pow(targetFP, 1.0 / Q);
    // std::cerr << "similarButUnequal = " << similarButUnequal << ", targetp1 = " << targetP1 << std::endl;
    // std::cerr << "UnsimilarButEqual = " << UnsimilarButEqual << ", targetp2 = " << targetP2 << std::endl;
    
    double l = 1e-6, r = gke.avgNGKmer * 2.0;
    double mid = (l + r) / 2.0;
    // while(l - r < - 1e-4) {
    //     mid = (r + l) / 2.0;
    //     if(__f(mid) - targetP1 < - 1e-4) {
    //         l = mid + 1e-6;
    //     } else {
    //         r = mid - 1e-6;
    //     }
    // }
    while(l - r < - 1e-4) {
        mid = (r + l) / 2.0;
        if(__rho(mid, Q1, Q80) - targetRho < - 1e-4) {
            r = mid + 1e-6;
        } else {
            l = mid - 1e-6;
        }
    }
    //return mid * 2.0 * gke.calcNGKmer(median_l); 
    mid *= gamma; 
    // return mid * gke.avgNGKmer * 2.0;
    return mid;
}


int PStableLSH::getOptimalQ(double sigma) {
    return std::ceil(std::log(node.n) / -std::log(__f(sigma / Q80 / gamma)));
}

int PStableLSH::getOptimalT(double sigma, int Q) {
    double p1_Q = 1.0, x = __f(sigma / Q1 / gamma);
    for(int i = 0; i < Q; i++) {
        p1_Q *= x;
    }
    return std::ceil(-std::log(targetFN[gke.k - 1])) / (1 - p1_Q);
}


double PStableLSH::getQ1(int w, int k, double l) {
    if(w == 1 && k == 1) {
        return -0.0004367 * l + 0.1575;
    } else if(w == 2 && k == 2) {
        return -0.000647 * l + 0.6234;
    } else if(w == 3 && k == 2) {
        return -0.001101 * l + 0.622;
    } else if(w == 3 && k == 3) {
        return 0.0001819 * l + 0.6861;
    } else if(w == 4 && k == 3) {
        return -7.222e-5 * l + 0.7186;
    } else if(w == 5 && k == 3) {
        return -0.0002641 * l + 0.7281;
    } else if(w == 5 && k == 4) {
        return 0.0002009 * l + 0.6702;
    }
    return 0.5;
}


double PStableLSH::getQ5(int w, int k, double l) {
    if(w == 1 && k == 1) {
        return -0.0006 * l + 0.2;
    } else if(w == 2 && k == 2) {
        return -0.0007714 * l + 0.6643;
    } else if(w == 3 && k == 2) {
        return -0.001203 * l + 0.6555;
    } else if(w == 3 && k == 3) {
        return 3.242e-5 * l + 0.7329;
    } else if(w == 4 && k == 3) {
        return -0.000204 * l + 0.7582;
    } else if(w == 5 && k == 3) {
        return -0.0003987 * l + 0.7656;
    } else if(w == 5 && k == 4) {
        return 2.203e-5 * l + 0.723;
    }
    return 0.5;
}

double PStableLSH::getQ20(int w, int k, double l) {
    if(w == 1 && k == 1) {
        return -0.00063 * l + 0.2225;
    } else if(w == 2 && k == 2) {
        return -0.000916 * l + 0.7102;
    } else if(w == 3 && k == 2) {
        return -0.001308 * l + 0.6916;
    } else if(w == 3 && k == 3) {
        return -8.313e-5 * l + 0.7747;
    } else if(w == 4 && k == 3) {
        return -0.0003313 * l + 0.7978;
    } else if(w == 5 && k == 3) {
        return -0.0005548 * l + 0.8093;
    } else if(w == 5 && k == 4) {
        return -0.0001776 * l + 0.7824;
    }
    return 0.5;
}


double PStableLSH::getQ80(int w, int k, double l) {
    if(w == 1 && k == 1) {
        return -0.00092 * l + 0.31;
    } else if(w == 2 && k == 2) {
        return -0.001056 * l + 0.7816;
    } else if(w == 3 && k == 2) {
        return -0.0015 * l + 0.7599;
    } else if(w == 3 && k == 3) {
        return -0.0004936 * l + 0.8893;
    } else if(w == 4 && k == 3) {
        return -0.0006091 * l + 0.8832;
    } else if(w == 5 && k == 3) {
        return -0.0008212 * l + 0.8874;
    } else if(w == 5 && k == 4) {
        return -0.0005698 * l + 0.8973;
    }
    return 0.5;
}


double PStableLSH::getQ95(int w, int k, double l) {
    if(w == 1 && k == 1) {
        return -0.00097 * l + 0.3375;
    } else if(w == 2 && k == 2) {
        return -0.001201 * l + 0.8275;
    } else if(w == 3 && k == 2) {
        return -0.00158 * l + 0.7922;
    } else if(w == 3 && k == 3) {
        return -0.000643 * l + 0.9361;
    } else if(w == 4 && k == 3) {
        return -0.0007019 * l + 0.9149;
    } else if(w == 5 && k == 3) {
        return -0.0008982 * l + 0.9151;
    } else if(w == 5 && k == 4) {
        return -0.0006897 * l + 0.9385;
    }
    return 0.5;
}


double PStableLSH::genCauchy() {
    std::uniform_real_distribution<double> uniDist(0.0, 1.0);
    std::random_device rd;
    // std::mt19937 gen(rd());
    std::mt19937 gen(nextRdSeed());
    return std::tan(std::acos(-1.0) * (uniDist(gen) - 0.5)) * gamma;
}


double PStableLSH::genUniInSigma() {
    std::uniform_real_distribution<double> uniDist(0.0, sigma);
    std::random_device rd;
    // std::mt19937 gen(rd());
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
            if(gke.k == 1) {
                HASHED_KMER gkmer = animoAcidHash(s[p]);
                if(a.count(gkmer) == 0) {
                    a[gkmer] = std::vector <double> (Q);
                    for(int k = 0; k < Q; ++ k) {
                        a[gkmer][k] = genCauchy();
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

        ///#pragma omp parallel for num_threads(16)
        for(int i = 0; i < node.n; ++ i) {
            std::vector <double> h(Q);
            for(int k = 0; k < Q; ++ k) {
                h[k] = 0.0;
            }

            auto& id = node.idList[i];
            auto& s = node.seqList.data[id];
            const int& l = s.length();

            std::map <HASHED_KMER, int> gkmerNumber;
            std::unordered_map <HASHED_KMER, std::string> gkmerToString;
            for(int p = 0; p < l; ++ p) {
                for(auto& mask: gke.otherIndices) {
                    if(p + mask.back() < l) {
                        HASHED_KMER gkmer = animoAcidHash(s[p]);
                        std::string kmer;
                        kmer += s[p];
                        for (auto& j: mask) {
                            gkmer <<= 5;
                            gkmer ^= animoAcidHash(s[p+j]);
                            kmer += s[p+j];
                        }
                        gkmerNumber[gkmer] ++;
                        gkmerToString[gkmer] = kmer;
                        /// for(int k = 0; k < Q; ++ k) {
                        ///     h[k] += a[gkmer][k];
                        /// }
                    }
                }
                if(gke.k == 1) {
                    HASHED_KMER gkmer = animoAcidHash(s[p]);
                    std::string kmer;
                    kmer += s[p];
                    gkmerNumber[gkmer] ++;
                    gkmerToString[gkmer] = kmer;
                }
            }

            for(int k = 0; k < Q; ++ k) {
                for(auto x: gkmerNumber) {
                    h[k] += a[x.first][k] * x.second;
                    // if(i < 3) {
                    //     std::cout << gkmerToString[x.first] << '\t' << a[x.first][k] << '\t' << x.second << std::endl;
                    // }
                }
                // if(i < 3) std::cout << std::endl;
            }
            
            // std::cout << i << ": ";
            for(int k = 0; k < Q; ++ k) {
                //h[k] /= gke.nGKmer[i] + gke.minNGKmer;
                //h[k] /= 2 * gke.nGKmer[i];
                // std::cout << h[k] << "\t";
                table[i].key[k] = (int)((h[k] + b[k]) / sigma);
            }
            // std::cout << std::endl;
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


    //    std::cout << "~~~ Glance DSU in" << T << std::endl;
    //    for(auto v: glanceDSU) {
    //        std::cout << v.first << ": ";
    //        for(auto x: v.second) {
    //            std::cout << x << ", ";
    //        }
    //        std::cout << std::endl;
    //    }
        
    }

    //Scan the dsu, and parse it to diffent sequence idlist
    subIdList = std::unordered_map <int, std::vector<int> > ();
    // subIdList = std::map <int, std::vector<int> > ();
    for(int i = 0; i < node.n; ++ i) {
        if(subIdList.count(dsu.find(i)) == 0) {
            subIdList[dsu.find(i)] = std::vector <int> ();
        }
        subIdList[dsu.find(i)].emplace_back(node.idList[i]);
    }
    if(subIdList.size() > 1) {
        for(auto x: subIdList) {
            for(auto y: x.second) {
                std::cout << y << " ";
            } 
            std::cout << std::endl;
        }
    }
}

