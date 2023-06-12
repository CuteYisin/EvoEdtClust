#include <cstring>
#include "HomoSearch.h"

void showUsage(std::string name) {
	std::cerr << "Usage: " << name << " <option(s)>\n"
		<< "Options:\n"
		<< "\t-h,\t\tShow this help message\n"
		<< "\t-i <input FASTA file path>\n"
		<< "\t-o <output ID pairs file path>\n"
		<< "\t-t,\t\tSimilarity threshold for two protein sequences to be considered homologous\n"
		<< std::endl;
}


int main(int argc, char* argv[]) {
    std::string inputFile = "", outputFile = "", threshold = "";
    double thresholdNum;
    if(argc < 6) {
		showUsage(argv[0]);
		return 1;
	}

    for (int i = 1; i < argc; ++ i) {
		if(strcmp(argv[i], "-h") == 0) {
			showUsage(argv[0]);
			return 0;
		} else if (strcmp(argv[i], "-o") == 0) {
			if (i + 1 < argc) {
				outputFile = argv[++i];
			} else {
				std::cerr << "-o option requires one argument." << std::endl;
				return 1;
			}
		} else if (strcmp(argv[i], "-i") == 0) {
			if (i + 1 < argc) {
				inputFile = argv[++i];
			} else {
				std::cerr << "-i option requires one argument." << std::endl;
				return 1;
			}
		} else if (strcmp(argv[i], "-t") == 0) {
			if (i + 1 < argc) {
				threshold = argv[++i];
                try {
                    thresholdNum = std::stod(threshold);
                } catch(const std::invalid_argument& e) {
                    std::cerr << "Error: -t invalid threshold." << std::endl;
                    return 1;
                }
			} else {
				std::cerr << "-t option requires one argument." << std::endl;
				return 1;
			}
		}
	}

	HomoSearch groundtruth;
	groundtruth.seqReservation(inputFile);
	groundtruth.filter(thresholdNum);
	groundtruth.printResult(outputFile);
    
    return 0;
}