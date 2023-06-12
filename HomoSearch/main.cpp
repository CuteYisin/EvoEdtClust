#include <cstring>
#include "Sequence.h"
#include "GappedKmerScan.h"
#include "HomoSearch.h"

void showUsage(std::string name) {
	std::cerr << "Usage: " << name << " <option(s)> input_fasta_file output\n"
		<< "Options:\n"
		<< "\t-h,\t\tShow this help message\n"
		<< "\t-c,\t\tList of parameters\n"
		<< std::endl;
}

int main(int argc, char* argv[]) {
	std::vector <int> windowSize, kmerLength, nRepetition;
	std::vector <double> projectionWidth, similarityThreshold;
	std::string inputFasta = "", output = "";

	if(argc < 5) {
		showUsage(argv[0]);
		return 1;
	}

	for (int i = 1; i < argc; ++ i) {
		if(strcmp(argv[i], "-h") == 0) {
			showUsage(argv[0]);
			return 0;
		} else if (strcmp(argv[i], "-c") == 0) {
			if (i + 1 < argc) {
				std::ifstream ifs(argv[++i]);
				if(ifs.is_open()) {
					while(!ifs.eof()) {
						int w, k, r;
						double p, s;
						ifs >> w >> k >> r >> p >> s;
						windowSize.emplace_back(w);
						kmerLength.emplace_back(k);
						nRepetition.emplace_back(r);
						projectionWidth.emplace_back(p); 
						similarityThreshold.emplace_back(s);
					}
				} else {
					std::cerr << "!!! Error: failed to open " << argv[i] << std::endl;
				}
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
        
	GappedKmerScan scanner(windowSize, kmerLength, nRepetition, projectionWidth);
	scanner.scan(seqList);

	LSH lsh(seqList, scanner, similarityThreshold);
	lsh.work();
	lsh.dumpToFile(output + "/output.pairs");

	return 0;
}
