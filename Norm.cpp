#include "Norm.h"


Norm::Norm(const GappedKmerScan& scanner): scanner(scanner) {}


Norm::~Norm() {}


L1::L1(const GappedKmerScan& scanner): Norm(scanner) {}


L1::~L1() {}


void L1::calculate(SequencePair& pairs) {
    int seqNumber = scanner.profile.size();
    for(int i = 0; i < seqNumber; i ++) {
        for(int j = i + 1; j < seqNumber; j ++) {
            double distance = 0;
            for(auto x: scanner.profile[i]) {
                if(scanner.profile[j].count(x.first) == 0) {
                    distance += x.second;
                } else {
                    distance += std::abs(x.second - scanner.profile[j].at(x.first));
                }
            }
            for(auto y: scanner.profile[j]) {
                if(scanner.profile[i].count(y.first) == 0) {
                    distance += y.second;
                }
            }
            pairs.addPair(i, j, distance, pairs.L1Distances);
        }
    }
}


L2::L2(const GappedKmerScan& scanner): Norm(scanner) {}


L2::~L2() {}


void L2::calculate(SequencePair& pairs) {
    int seqNumber = scanner.profile.size();
    for(int i = 0; i < seqNumber; i ++) {
        for(int j = i + 1; j < seqNumber; j ++) {
            double distance = 0;
            for(auto x: scanner.profile[i]) {
                if(scanner.profile[j].count(x.first) == 0) {
                    distance += x.second * x.second;
                } else {
                    double value = std::abs(x.second - scanner.profile[j].at(x.first));
                    distance += value * value;
                }
            }
            for(auto y: scanner.profile[j]) {
                if(scanner.profile[i].count(y.first) == 0) {
                    distance += y.second * y.second;
                }
            }
            distance = sqrt(distance);
            pairs.addPair(i, j, distance, pairs.L2Distances);
        }
    }
}


Cosine::Cosine(const GappedKmerScan& scanner): Norm(scanner) {}


Cosine::~Cosine() {}


void Cosine::calculate(SequencePair& pairs) {
    int seqNumber = scanner.profile.size();
    for(int i = 0; i < seqNumber; i ++) {
        for(int j = i + 1; j < seqNumber; j ++) {
            double dotProduct = 0, normA = 0, normB = 0, distance = 0;
            for(auto x: scanner.profile[i]) {
                normA += x.second * x.second;
                if(scanner.profile[j].count(x.first)) {
                    dotProduct += x.second * scanner.profile[j].at(x.first);
                }
            }
            for(auto y: scanner.profile[j]) {
                normB += y.second * y.second;
            }
            distance = 1.0 - dotProduct / sqrt(normA * normB);
            pairs.addPair(i, j, distance, pairs.cosineDistances);
        }
    }
}


Jaccard::Jaccard(const GappedKmerScan& scanner): Norm(scanner) {}


Jaccard::~Jaccard() {}


void Jaccard::calculate(SequencePair& pairs) {
    int seqNumber = scanner.profile.size();
    for(int i = 0; i < seqNumber; i ++) {
        for(int j = i + 1; j < seqNumber; j ++) {
            double intersection = 0, sum = 0, distance = 0;
            for(auto x: scanner.profile[i]) {
                if(scanner.profile[j].count(x.first)) {
                    intersection += std::min(x.second, scanner.profile[j].at(x.first));
                    sum += x.second;
                }
            }
            for(auto y: scanner.profile[j]) {
                sum += y.second;
            }
            distance = 1.0 - intersection / (sum - intersection);
            pairs.addPair(i, j, distance, pairs.jaccardDistances);
        }
    }
}