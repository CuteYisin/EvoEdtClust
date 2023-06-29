#include "Sequence.h"


SequenceList::SequenceList() {
    std::cout << "\n=== Collecting a list of sequences " << std::endl;
    data = std::vector <std::string> ();
}


SequenceList::~SequenceList() { }


void SequenceList::loadFromFasta_dumpIndexedHeader(const std::string& fastaFile, const std::string& outputFile) {
    std::cout << "+++ Load fasta file from " << fastaFile << std::endl;
    std::ifstream ifs(fastaFile);
    if(!ifs.is_open()) {
        std::cerr << "!!! Error: failed to open " << fastaFile << std::endl;
    }

    std::cout << "+++ Dump header indices to file: " << outputFile << std::endl;
    std::ofstream ofs(outputFile);
    if(!ofs.is_open()) {
        std::cerr << "!!! Error: failed to open " << outputFile << std::endl;
    }

    std::string header, sequence;
    int index = 0;
    while (!ifs.eof()) {
        std::getline(ifs, header);
        if (header[0] == '>') {
            std::getline(ifs, sequence);
            data.emplace_back(sequence);
            //ofs << index << "\t" << header << "\n";
            ofs << header << "\n";
            ++ index;
        }
    }
    ofs.flush();
    ofs.close();
    std::cout << "=== Finally, " << index << " sequences loaded." << std::endl;
}


ClusterNode::ClusterNode(const SequenceList& seqList, 
                         const std::vector <int>& idList, 
                         int level, double expectedSimLowBound):
    seqList(seqList), level(level), expectedSimLowBound(expectedSimLowBound) { 
        this->idList = std::vector <int> ();
        n = 0, avgL = 0.0;
        for(auto& id: idList) {
            this->idList.emplace_back(id);
            ++ n, avgL += seqList.data[id].length();
        }
        avgL /= n;
}


ClusterNode::~ClusterNode() { }


void ClusterNode::show() {
    std::cout << "=== Current set contains " << n << " sequences with average length " << avgL << " .. " << std::endl; 
    std::cout << "+++ Sequences are: ";
    int glanceN = 10;
    if((int)idList.size() <= glanceN) {
        for(auto id: idList) {
            std::cout << id << ", ";
        }
        std::cout << "." << std::endl;
    } else {
        for(int i = 0; i < glanceN/2; i ++) {
            std::cout << idList[i] << ", ";
        }
        std::cout << "..." ;
        for(int i = 0; i < glanceN/2; i ++) {
            std::cout << ", " << idList[idList.size() - glanceN/2 + i];
        }
        std::cout << "." << std::endl;
    }
    std::cout << "+++ This set is level-" << level << " partition" << std::endl; 
    std::cout << "+++ The expected similarity lower bound in this set is " << expectedSimLowBound << std::endl; 
}


ClusterTree::ClusterTree(const SequenceList& seqList) {
    std::vector <int> idList;
    for(int i = 0; i < (int)seqList.data.size(); ++ i) {
        idList.emplace_back(i);
    }

    bfsOrder = std::queue <ClusterNode> ();
    bfsOrder.push(ClusterNode(seqList, idList, 0, 0.0));
}


ClusterTree::~ClusterTree() { }