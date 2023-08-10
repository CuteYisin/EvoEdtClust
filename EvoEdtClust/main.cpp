#include <cstring>
#include <ctime>

#include "Sequence.h"
#include "GappedKmer.h"
#include "LSH.h"


void showUsage(std::string name) {
	std::cerr << "Usage: " << name << " <option(s)> input_fasta_file output\n"
		<< "Options:\n"
		<< "\t-h,\t\tShow this help message\n"
		<< "\t-t,\t\tthreshold of similarity (default: 0.5)\n"
		<< std::endl;
}


int main(int argc, char* argv[]) {
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(NULL);
	std::cout.tie(NULL);
	std::cerr.tie(NULL);

	std::string inputSim = "";
	double similarityThreshold = 0.5;
	std::string inputFastaFile = "", outputDir = "";

	if(argc < 2) {
		showUsage(argv[0]);
		return 1;
	}

	for (int i = 1; i < argc; ++ i) {
		if(strcmp(argv[i], "-h") == 0) {
			showUsage(argv[0]);
			return 0;
		} else if (strcmp(argv[i], "-t") == 0) {
			if (i + 1 < argc) {
				inputSim = argv[++i];
				similarityThreshold = std::stof(inputSim);
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
	
    SequenceList seqList;
	std::string headerIndicesFile = outputDir + "/headerIndices.out";
	seqList.loadFromFasta_dumpIndexedHeader(inputFastaFile, headerIndicesFile);
        
	ClusterTree tree(seqList);
	std::vector <ClusterNode> leaves;

	while(!tree.bfsOrder.empty()) {
		auto& node = tree.bfsOrder.front();
		if(node.n > 10 && node.level< (int)ParameterGenerator::w.size()) {
			int w = ParameterGenerator::w[node.level];
			int k = ParameterGenerator::k[node.level];
			double sim = ParameterGenerator::computeSimFromWK(w, k, node.avgL);

			std::cout << "+++ w = " << w << ", k = " << k << ", sim = " << sim << std::endl;

			GappedKmerEmbedding gkmer(node, w, k);
			gkmer.scan();

			PStableLSH plsh(node, gkmer);
			plsh.scanPars();
			plsh.work();
			std::cout << std::endl;
			for(auto& x: plsh.subIdList) {
				auto& idList = x.second;
				ClusterNode newNode(seqList, idList, node.level + 1, sim);
				tree.bfsOrder.push(newNode);
			}
		} else {
			leaves.push_back(node);
		}
		tree.bfsOrder.pop();
	}

	std::cout << "===\n===\n=== " << std::endl;
	long long nPairs = 0;
	for(auto& leaf: leaves) {
		nPairs += leaf.n * (leaf.n-1) / 2;
	}

	std::cout << "=== Finally, " << leaves.size() << " clusters are found.. " << std::endl;
	std::cout << "=== Finally, " << nPairs << " pairs are found.. " << std::endl;

	
	std::vector <std::string> headerList;
	std::ifstream ifs(headerIndicesFile);
	std::string header;
	while(!ifs.eof()){
		std::getline(ifs, header);
		headerList.emplace_back(header);
	}

	std::string finalResFile = outputDir + "/EvoEdtClust_" + inputSim + ".out";

    std::cout << "+++ Dump clusters to file: " << finalResFile << std::endl;
    std::ofstream ofs(finalResFile);
    if(!ofs.is_open()) {
        std::cerr << "!!! Error: failed to open " << finalResFile << std::endl;
    }

	int nCluster = 0;
	for(auto& leaf: leaves) {
		if(leaf.n >= 2) {
			nCluster ++;
			ofs << "#_Cluster_" << nCluster << "\n";
			for(int i = 0; i < leaf.n; i ++) {
				auto& id = leaf.idList[i];
				// ofs << seqList.data[id] << "\n";
				ofs << headerList[id] << "\n";
			}
		}
	}
	ofs.flush();
	ofs.close();

	return 0;
}
