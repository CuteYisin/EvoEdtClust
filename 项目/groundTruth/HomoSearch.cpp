#include "HomoSearch.h"

HomoSearch::HomoSearch() {
    seq = std::vector<std::string> ();
    HomoPair = std::vector <std::pair <ID_TYPE, ID_TYPE>> ();
}


HomoSearch::~HomoSearch() { }


void HomoSearch::seqReservation(const std::string& inputFasta) {
    std::cout << "+++ Load Fasta file from " << inputFasta << std::endl;

    std::ifstream ifs(inputFasta);
    std::string header, sequence;
    ID_TYPE progress = 0;
    if(ifs.is_open()) {
        while(std::getline(ifs, header) && std::getline(ifs, sequence)) {
            seq.emplace_back(sequence);
            progress ++;
        }
        std::cout << "=== Finally, " << progress << " sequences loaded." << std::endl;
    } else {
        std::cerr << "!!! Error: failed to open " << inputFasta << std::endl;
    }
}


double HomoSearch::EditSimilarity(const std::string& str1, const std::string& str2) {
    int len1 = str1.size(), len2 = str2.size();
    if(! len1 || ! len2) {
        return 0.0;
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

    return 1 - (double)dp[len1][len2] / std::max(len1, len2);
}


void HomoSearch::filter(const double& threshold) {
    std::cout << "=== Start generating the ground truth..." << std::endl;
    for(std::vector<std::string>::size_type i = 0; i < seq.size(); ++i) {
        for(std::vector<std::string>::size_type j = i + 1; j < seq.size(); ++j) {
                if(EditSimilarity(seq[i], seq[j]) - threshold > 1e-8) {
                    HomoPair.emplace_back(std::make_pair(i, j));
            }
        }
    }
}


void HomoSearch::printResult(const std::string& outputFile) {
    std::ofstream ofs(outputFile);
	if(ofs.is_open()) {
		for(const auto& Pair : HomoPair) {
            ofs << std::min(Pair.first, Pair.second) << '\t' << std::max(Pair.first, Pair.second) << std::endl;
        }
		ofs.close();
	} else {
		std::cerr << "!!! Error: failed to open " << outputFile << std::endl;
	}
}
