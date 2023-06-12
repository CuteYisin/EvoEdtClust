#include "EditDistance.h"

EditDistance::EditDistance(const SequenceList& seqList) :seqList(seqList) {}


EditDistance::~EditDistance() {}


int EditDistance::getEditDistance(const std::string& str1, const std::string& str2) {
    int len1 = str1.size(), len2 = str2.size();
    if(! len1 || ! len2) {
        return len1 + len2;
    }

    std::vector<std::vector<int>> dp(len1+1, std::vector<int>(len2+1, 0));
    for(int i = 0; i <= len1; ++i) {
        dp[i][0] = i;
    }
    for(int j = 0; j <= len2; ++j) {
        dp[0][j] = j;
    }

    for(int i = 1; i <= len1; ++i) {
        for(int j = 1; j <= len2; ++j) {
            if(str1[i-1] == str2[j-1]) {
                dp[i][j] = dp[i-1][j-1];
            } else {
                dp[i][j] = std::min(dp[i-1][j-1], std::min(dp[i-1][j], dp[i][j-1])) + 1;
            }
        }
    }

    return dp[len1][len2];
}


void EditDistance::calculate(SequencePair& pairs) {
    int seqNumber = seqList.list.size();
    double distance = 0;
    for(int i = 0; i < seqNumber; i ++) {
        for(int j = i + 1; j < seqNumber; j ++) {
            distance = getEditDistance(seqList.list[i].sequence, seqList.list[j].sequence);
            pairs.addPair(i, j, distance, pairs.editDistances);
        }
    }
}