#include <cstring>
#include <ctime>
#include "Sequence.h"
#include "GappedKmer.h"
#include "LSH.h"


void showUsage(std::string name) {
	std::cerr << "Usage: " << name << " <option(s)> input_fasta_file output\n"
		<< "Options:\n"
		<< "\t-h,\t\tShow this help message\n"
		<< "\t-m,\t\tSelect one of LSH Methods (1. 1-Stable distributions (default); 2. MinHash)\n"
		<< "\t-t,\t\tthreshold of similarity (default: 0.5)\n"
		<< std::endl;
}

int main(int argc, char* argv[]) {
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(NULL);
	std::cout.tie(NULL);
	std::cerr.tie(NULL);

	int mode = 1;
	double similarityThreshold = 0.5;
	std::string inputFastaFile = "", outputDir = "";

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
				std::cerr << "-m option requires one argument." << std::endl;
				return 1;
			}  
		} else if (strcmp(argv[i], "-t") == 0) {
			if (i + 1 < argc) {
				similarityThreshold = std::stof(argv[++i]);
			} else {
				std::cerr << "-t option requires one argument." << std::endl;
				return 1;
			}  
		} else {
			if(i + 1 < argc) {
				inputFastaFile = argv[i++];
				outputDir = argv[i];
			} else {
				std::cerr << "require input file path and output file directory";
			}
		}
	}
	
	//auto st = clock();
    //SequenceList seqList;
    //seqList.loadFromFasta(inputFastaFile);
	//seqList.dumpIndexedHeader(outputDir + "/headerIndice.out");
	//auto ed = clock();
	//std::cerr << "~~~ Load data time usage " << (double)(ed - st ) / CLOCKS_PER_SEC << std::endl;
	auto st = clock();
    LiteSequenceList seqList;
	seqList.loadFromFasta_dumpIndexedHeader(inputFastaFile, outputDir + "/headerIndices.out");
	auto ed = clock();
	std::cerr << "~~~ Load data time usage " << (double)(ed - st ) / CLOCKS_PER_SEC << std::endl;
        
	/*
	for(int i = 0; i < (int)seqList.list.size(); i ++) {
		std::cout << seqList.list[i] << std::endl;
	}
	*/
	
	st = clock();
	SequenceClusterTree seqClusterTree(seqList);
	std::vector <SequenceClusterNode> seqClusterTreeLeaves;
	while(!seqClusterTree.treeInBFS.empty()) {
		SequenceClusterNode& seqClusterNode = seqClusterTree.treeInBFS.front();
		if(seqClusterNode.n >= 10 && seqClusterNode.level <= 0 &&
			seqClusterNode.expectedSimLowBound - similarityThreshold < -1e-6) {
			//std::cout << std::endl;
			//seqClusterNode.show();

			double nextSimLowBound = seqClusterNode.expectedSimLowBound;
			int w = 0, k = 0;
			if(seqClusterNode.level == 0) {
				w = 2, k = 2;
				nextSimLowBound = ParameterGenerator::computeSimFromWK(w, k, seqClusterNode.avgL);
			} else {
				nextSimLowBound = std::min(similarityThreshold, seqClusterNode.expectedSimLowBound + 0.1);
				ParameterGenerator::updateWKFromSim(nextSimLowBound, w, k, seqClusterNode.avgL);
			}

			std::cout << "+++ Now wish partition this set in to cluster with expected similariy " << nextSimLowBound << std::endl;
			std::cout << "+++ w = " << w << ", k = " << k << std::endl;

			auto s1 = clock();
			GappedKmerEmbedding gkmerEmbedding(w, k);
			gkmerEmbedding.scan(seqClusterNode);
			std::cerr << "~~~ GappedKmerEmbedding data time usage " << (double)(clock() - s1) / CLOCKS_PER_SEC << std::endl;

			s1 = clock();
			PStableLSH plsh(seqClusterNode, gkmerEmbedding);
			plsh.work();
			std::cerr << "~~~ GappedKmerEmbedding data time usage " << (double)(clock() - s1) / CLOCKS_PER_SEC << std::endl;

			for(auto& sub: plsh.finalSubset) {
				SequenceClusterNode newCluster(seqList, sub.second, seqClusterNode.level + 1, nextSimLowBound);
				//newCluster.show();
				seqClusterTree.treeInBFS.push(newCluster);
			}
		} else {
			seqClusterTreeLeaves.push_back(seqClusterNode);
		}
		seqClusterTree.treeInBFS.pop();
	}
	ed = clock();
	std::cerr << "~~~ BuildCluster data time usage " << (double)(ed - st ) / CLOCKS_PER_SEC << std::endl;

	std::cout << "=== === === " << std::endl;
	for(auto& seqClusterLeaf: seqClusterTreeLeaves) {
		seqClusterLeaf.show();
	}
	std::cout << "=== Finally, " << seqClusterTreeLeaves.size() << " clusters are found.. " << std::endl;

	return 0;
}
