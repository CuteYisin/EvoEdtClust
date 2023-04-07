#include <iostream>
#include <string>
#include <memory>
#include "SeqCode.h"
#include "SparseMatrixMaker.h"

void showUsage(std::string name) {
	std::cerr << "Usage: " << name << " <option(s)>\n"
		<< "Options:\n"
		<< "\t-h,\t\tShow this help message\n"
		<< "\t-m <Linclust|ALFATClust|LSH>,\t\tThe mode based on different algorithm inputs\n"
		<< "\t-l <LookupTable file path>\n"
		<< "\t-i <cluster file path>\n"
		<< "\t-o <output sparse-matrix file path>\n"
		<< std::endl;
}


int main(int argc, char* argv[]){
    std::string mode = "", lookupTable = "", inputFile = "", outputFile = "";
    if(argc < 8) {
		showUsage(argv[0]);
		return 1;
	}

    for (int i = 1; i < argc; ++ i) {
		if(strcmp(argv[i], "-h") == 0) {
			showUsage(argv[0]);
			return 0;
		} else if (strcmp(argv[i], "-m") == 0) {
			if (i + 1 < argc) {
				mode = argv[++i];
			} else {
				std::cerr << "-m option requires one argument." << std::endl;
				return 1;
			}
		} else if (strcmp(argv[i], "-l") == 0) {
			if (i + 1 < argc) {
				lookupTable = argv[++i];
			} else {
				std::cerr << "-l option requires one argument." << std::endl;
				return 1;
			}
		} else if (strcmp(argv[i], "-i") == 0) {
			if (i + 1 < argc) {
				inputFile = argv[++i];
			} else {
				std::cerr << "-i option requires one argument." << std::endl;
				return 1;
			}
		} else if (strcmp(argv[i], "-o") == 0) {
			if (i + 1 < argc) {
				outputFile = argv[++i];
			} else {
				std::cerr << "-o option requires one argument." << std::endl;
				return 1;
			}
		}
	}

    SeqCode table;
    table.code(lookupTable);

    std::unique_ptr<SparseMatrixMaker> Conv;
    if(mode == "Linclust") {
        Conv = std::make_unique<LinclustConv>(table);
    } else if(mode == "ALFATClust") {
        Conv = std::make_unique<ALFATConv>(table);
    } else if(mode == "LSH") {
        Conv = std::make_unique<LSHConv>(table);
    } else {
        std::cerr << "-m option is provided wrong argument." << std::endl;
        return 1;
    }

    Conv->convert(inputFile);
    Conv->outputAsAllPairs(outputFile);

    return 0;
}
