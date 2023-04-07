#include "SparseMatrixMaker.h"


HomoPair::HomoPair() {
	seqPair = std::unordered_map <std::string, double> ();
}


HomoPair::~HomoPair() { }


inline std::string HomoPair::indexPairToString(const ID_TYPE& i, const ID_TYPE& j) {
	return std::to_string(i) + "-" + std::to_string(j);
}


inline std::pair <ID_TYPE, ID_TYPE> HomoPair::stringToIndexPair(const std::string& s) {
    int pos = s.find('-');
    return std::make_pair(std::stoi(s.substr(0, pos)), std::stoi(s.substr(1+pos, s.size()-pos-1)));
}


void HomoPair::create(const ID_TYPE& i, const ID_TYPE& j, const double& sim) {
    seqPair[indexPairToString(i, j)] = sim;
}


SparseMatrixMaker::SparseMatrixMaker(const SeqCode& table): table(table), cluster(std::vector<ID_TYPE>()) {}


SparseMatrixMaker::~SparseMatrixMaker() {}


void SparseMatrixMaker::convert(const std::string& clusterFile) {}


void SparseMatrixMaker::outputAsAllPairs(const std::string& outputFilePath) {
    std::ofstream ofs(outputFilePath);
	if(ofs.is_open()) {
        for (const auto& elem : result.seqPair) {
            if(elem.second) {
                auto tmp = HomoPair::stringToIndexPair(elem.first);
                ofs << tmp.first << '\t' << tmp.second << '\t' << elem.second <<std::endl;
            }
        }
		ofs.close();
	} else {
		std::cerr << "!!! Error: failed to open " << outputFilePath << std::endl;
	}
}


LinclustConv::LinclustConv(const SeqCode& table): SparseMatrixMaker(table) {
    std::cout << "--- Start making sparse matrix of Linclust" << std::endl;
}


LinclustConv::~LinclustConv() {}


void LinclustConv::convert(const std::string& clusterFile) {
    std::cout << "+++ Load cluster file from " << clusterFile << std::endl;

    std::ifstream ifs(clusterFile);
    std::string tmp, header, sequence;
    int progress = 0;
    if(ifs.is_open()) {
		while(!ifs.eof()) {
            std::vector <ID_TYPE> ().swap(cluster);

            if(tmp == "") {
                std::getline(ifs, tmp);
            } else {
                header = tmp;
                std::getline(ifs, sequence);
                cluster.emplace_back(table.seqID.at(header));
            }
            while(std::getline(ifs, header) && std::getline(ifs, sequence)) {
                if(sequence[0] == '>') {
                    tmp = sequence;
                    break;
                }
                cluster.emplace_back(table.seqID.at(header));
            }

            for(auto i = cluster.begin(); i!= cluster.end(); i++) {
                for(auto j = i+1; j!= cluster.end(); j++) {
                    result.create(std::min(*i, *j), std::max(*i, *j), 1);
                    progress ++;
                }
            }
		}
		std::cout << "=== Finally, the sparse matrix has " << progress << " non-empty elements." << std::endl;

	} else {
		std::cerr << "!!! Error: failed to open " << clusterFile << std::endl;
	}
}


ALFATConv::ALFATConv(const SeqCode& table): SparseMatrixMaker(table) {
    std::cout << "--- Start making sparse matrix of ALFATClust" << std::endl;
}


ALFATConv::~ALFATConv() {}


void ALFATConv::convert(const std::string& clusterFile) {
    std::cout << "+++ Load cluster file from " << clusterFile << std::endl;

    std::ifstream ifs(clusterFile);
    std::string tmp, header;
    int progress = 0;
    if(ifs.is_open()) {
		while(!ifs.eof()) {
            std::vector <ID_TYPE> ().swap(cluster);

            if(tmp == "") {
                std::getline(ifs, tmp);
            }
            while(std::getline(ifs, header)) {
                if(header[0] == '#') {
                    tmp = header;
                    break;
                }
                header.insert(header.begin(), '>');
                cluster.emplace_back(table.seqID.at(header));
            }

            for(auto i = cluster.begin(); i!= cluster.end(); i++) {
                for(auto j = i+1; j!= cluster.end(); j++) {
                    result.create(std::min(*i, *j), std::max(*i, *j), 1);
                    progress ++;
                }
            }
		}
		std::cout << "=== Finally, the sparse matrix has " << progress << " non-empty elements." << std::endl;

	} else {
		std::cerr << "!!! Error: failed to open " << clusterFile << std::endl;
	}
}


LSHConv::LSHConv(const SeqCode& table): SparseMatrixMaker(table) {
    std::cout << "--- Start making sparse matrix of LSH" << std::endl;
}


LSHConv::~LSHConv() {}


void LSHConv::convert(const std::string& clusterFile) {
    std::cout << "+++ Load cluster file from " << clusterFile << std::endl;

    std::ifstream ifs(clusterFile);
    std::string headerA, headerB;
    double sim;
    int progress = 0;
    if(ifs.is_open()) {
		while(!ifs.eof()) {
            std::getline(ifs, headerA, '\t');
            std::getline(ifs, headerB, '\t');
            ifs >> sim;
            ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            int IDA = table.seqID.at(headerA);
            int IDB = table.seqID.at(headerB);

            result.create(std::min(IDA, IDB), std::max(IDA, IDB), sim);
            progress ++;
		}
		std::cout << "=== Finally, the sparse matrix has " << progress << " non-empty elements." << std::endl;

	} else {
		std::cerr << "!!! Error: failed to open " << clusterFile << std::endl;
	}
}
