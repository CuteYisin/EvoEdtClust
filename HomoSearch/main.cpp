#include <cstring>
#include "Sequence.h"
#include "GappedKmerScan.h"
#include "HomoSearch.h"

void showUsage(std::string name) {
	std::cerr << "Usage: " << name << " <option(s)> input_fasta_file output\n"
		<< "Options:\n"
		<< "\t-h,\t\tShow this help message\n"
		<< "\t-m,\t\tSelect one of LSH Methods (1. 1-Stable distributions (default); 2. MinHash)\n"
		<< "\t-t,\t\tthreshold of similarity (default: 0.5)\n"
		<< std::endl;
}

int main(int argc, char* argv[]) {
	std::vector <int> windowSize = {2, 10, 5, 20, 30};
	std::vector <int> kmerLength = {2, 2, 3, 3, 3};
	int mode = 1;
	double similarityThreshold = 0.5;
	std::string inputFasta = "", output = "";

	if(argc < 3) {
		showUsage(argv[0]);
		return 1;
	}

	for (int i = 1; i < argc; ++ i) {
		if(strcmp(argv[i], "-h") == 0) {
			showUsage(argv[0]);
			return 0;
		} else if (strcmp(argv[i], "-m") == 0) {
			if (i + 1 < argc) {
				mode = std::stoi(argv[++i]);
			} else {
				std::cerr << "-c option requires one argument." << std::endl;
				return 1;
			}  
		} else if (strcmp(argv[i], "-t") == 0) {
			if (i + 1 < argc) {
				similarityThreshold = std::stof(argv[++i]);
			} else {
				std::cerr << "-c option requires one argument." << std::endl;
				return 1;
			}  
		} else {
			if(i + 1 < argc) {
				inputFasta = argv[i++];
				output = argv[i];
			} else {
				std::cerr << "require both input file path and output file path";
			}
		}
	}
	
    SequenceList seqList;
    seqList.loadFromFasta(inputFasta);
        
	GappedKmerScan scanner(windowSize, kmerLength);
	scanner.scan(seqList);

	LSH lsh(seqList, scanner, mode, similarityThreshold);
	lsh.work();
	lsh.dumpToFile(output + "/output.pairs");

	return 0;
}
