#include <cstring>
#include "Sequence.h"
#include "GappedKmerScan.h"
#include "Norm.h"
#include "EditDistance.h"
#include "LSH.h"

void showUsage(std::string name) {
	std::cerr << "Usage: " << name << " <option(s)> input_fasta_file output\n"
		<< "Options:\n"
		<< "\t-h,\t\tShow this help message\n"
		<< "\t-w,\t\tThe window size (default: 5)\n"
		<< "\t-k,\t\tThe kmer length (default: 3)\n"
		<< "\t-p,\t\tThe number of repeated hashing (default: 1)\n"
		<< std::endl;
}

int main(int argc, char* argv[]) {
	int windowSize = 5, kmerLength = 3, repetition = 1;
	std::string inputFasta = "", output = "", parameter="";

	if(argc < 3) {
		showUsage(argv[0]);
		return 1;
	}

	for (int i = 1; i < argc; ++ i) {
		if(strcmp(argv[i], "-h") == 0) {
			showUsage(argv[0]);
			return 0;
		} else if (strcmp(argv[i], "-w") == 0) {
			if (i + 1 < argc) {
				windowSize = std::stoi(argv[++i]); 
			} else {
				std::cerr << "-w option requires one argument." << std::endl;
				return 1;
			}  
		} else if (strcmp(argv[i], "-k") == 0) {
			if (i + 1 < argc) {
				kmerLength = std::stoi(argv[++i]); 
			} else {
				std::cerr << "-k option requires one argument." << std::endl;
				return 1;
			}  
		} else if (strcmp(argv[i], "-p") == 0) {
			if (i + 1 < argc) {
				repetition = std::stoi(argv[++i]); 
			} else {
				std::cerr << "-p option requires one argument." << std::endl;
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
	parameter = "/" + std::to_string(windowSize) + "-" + std::to_string(kmerLength);
	
    SequenceList seqList;
    seqList.loadFromFasta(inputFasta);

	SequencePair pairs(seqList);

	/*EditDistance Ed(seqList);
	Ed.calculate(pairs);
	pairs.dumpToFile(output + "/EditDistance.txt", pairs.editDistances);*/
        
	GappedKmerScan scanner(windowSize, kmerLength);
	scanner.scan(seqList);

	/*L1 l1(scanner);
	l1.calculate(pairs);
	pairs.dumpToFile(output + parameter + "/L1Distance.txt", pairs.L1Distances);

	L2 l2(scanner);
	l2.calculate(pairs);
	pairs.dumpToFile(output + parameter + "/L2Distance.txt", pairs.L2Distances);

	Cosine Cos(scanner);
	Cos.calculate(pairs);
	pairs.dumpToFile(output + parameter + "/CosineDistance.txt", pairs.cosineDistances);

	Jaccard jaccard(scanner);
	jaccard.calculate(pairs);
	pairs.dumpToFile(output + parameter + "/JaccardDistance.txt", pairs.jaccardDistances);*/

	pStable p(repetition, scanner);
	p.work(1);
	p.print(output + parameter + "/L1DistanceEstimation-" + std::to_string(repetition) + ".txt", 1);

	//p.work(2);
	//p.print(output + parameter + "/L2DistanceEstimation-" + std::to_string(repetition) + ".txt", 1);

	RandomProjection random(repetition, scanner);
	random.work();
	random.print(output + parameter + "/CosineDistanceEstimation-" + std::to_string(repetition) + ".txt", 2);

	MinHash minhash(repetition, scanner);
	minhash.work();
	minhash.print(output + parameter + "/JaccardDistanceEstimation-" + std::to_string(repetition) + ".txt", 3);

	return 0;
}
