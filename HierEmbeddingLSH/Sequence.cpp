#include "Sequence.h"


Sequence::Sequence(const std::string& header, const std::string& sequence, const int& index):
    header(header), sequence(sequence), index(index) { }


Sequence::~Sequence() { }


SequenceList::SequenceList() {
    std::cout << "\n=== Collecting a list of sequences " << std::endl;
    list = std::vector <Sequence> ();
}


SequenceList::~SequenceList() { }


void SequenceList::loadFromFasta(const std::string& fastaFile) {
    std::cout << "+++ Load fasta file from " << fastaFile << std::endl;

    std::ifstream ifs(fastaFile);
    std::string header, sequence; //can be optimized using string_view, which is an ultrafast read-only string type supported in C++17
    int index = 0;
    if(ifs.is_open()) {
        while(!ifs.eof()) {
            std::getline(ifs, header);
            if(header[0] == '>') {
                std::getline(ifs, sequence);
                list.emplace_back(Sequence(header, sequence, index));

                index ++;
                if(index % 10000 == 0) {
                    std::cout << "--- " << index << " sequences loaded ..." << std::endl;
                }
            }
        }
        std::cout << "+++ Finally, " << index << " sequences loaded." << std::endl;

    } else {
        std::cerr << "!!! Error: failed to open " << fastaFile << std::endl;
    }
}


void SequenceList::dumpIndexedHeader(const std::string& outputFile) {
    std::cout << "+++ Dump header indices to file: " << outputFile << std::endl;

    std::ofstream ofs(outputFile);
    if(ofs.is_open()) {
        for(auto s: list) {
            ofs << s.index << "\t" << s.header << std::endl;
        }
        ofs.flush();
        ofs.close();
    } else {
        std::cerr << "!!! Error: failed to open " << outputFile << std::endl;
    }
}


LiteSequenceList::LiteSequenceList() {
    std::cout << "\n=== Collecting a list of sequences " << std::endl;
    list = std::vector <std::string> ();
}


LiteSequenceList::~LiteSequenceList() { }


void LiteSequenceList::loadFromFasta_dumpIndexedHeader(const std::string& fastaFile, const std::string& outputFile) {
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

    /*
    std::string_view header, sequence;
    int index = 0;

    int fd = open(fastaFile.c_str(), O_RDONLY);
    struct stat statbuf;
    fstat(fd, &statbuf);
    const char *ptr = static_cast<const char *>(mmap(NULL, statbuf.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
    size_t res = statbuf.st_size;
    const char *p = ptr;
    const char *q = NULL;
    unsigned long length;
    int n = std::count(ptr, ptr + statbuf.st_size, '\n') / 2;
    while(index < n) {
        q = static_cast<const char *>(std::memchr(p, '\n', res));
        length = static_cast<int>(q - p);
        header = std::string_view({p, length});
        p = q + 1, res -= length;

        q = static_cast<const char *>(std::memchr(p, '\n', res));
        length = static_cast<int>(q - p);
        sequence = std::string_view({p, length});
        p = q + 1, res -= length;

        list.emplace_back(sequence);
        //ofs << index << "\t" << header << "\n";
        ++ index;
    }
    ofs.flush();
    ofs.close();
    */

    std::string header, sequence;
    int index = 0;
    while (!ifs.eof()) {
        std::getline(ifs, header);
        if (header[0] == '>') {
            std::getline(ifs, sequence);
            list.emplace_back(sequence);
            ofs << index << "\t" << header << "\n";
            ++index;
        }
    }
    ofs.flush();
    ofs.close();
    std::cout << "+++ Finally, " << index << " sequences loaded." << std::endl;
}


SequenceClusterNode::SequenceClusterNode(const LiteSequenceList& seqList, 
                                         const std::vector <int>& seqSet, 
                                         int level, double expectedSimLowBound):
    seqList(seqList), level(level), expectedSimLowBound(expectedSimLowBound) { 
        this->seqSet = std::vector <int> ();
        n = 0, avgL = 0.0;
        for(const auto& s: seqSet) {
            this->seqSet.emplace_back(s);
            ++ n, avgL += seqList.list[s].length();
        }
        avgL /= n;
}


SequenceClusterNode::~SequenceClusterNode() { }


void SequenceClusterNode::show() {
    std::cout << "=== Current set contains " << n << " sequences with average length " << avgL << " .. " << std::endl; 
    std::cout << "+++ Sequences are: ";
    int glanceN = 10;
    if((int)seqSet.size() <= glanceN) {
        for(auto i: seqSet) {
            std::cout << i << ", ";
        }
        std::cout << "." << std::endl;
    } else {
        for(int i = 0; i < glanceN/2; i ++) {
            std::cout << seqSet[i] << ", ";
        }
        std::cout << "..." ;
        for(int i = 0; i < glanceN/2; i ++) {
            std::cout << ", " << seqSet[seqSet.size() - glanceN/2 + i];
        }
        std::cout << "." << std::endl;
    }
    std::cout << "+++ This set is level-" << level << " partition" << std::endl; 
    std::cout << "+++ The expected similarity lower bound in this set is " << expectedSimLowBound << std::endl; 
}


SequenceClusterTree::SequenceClusterTree(const LiteSequenceList& seqList) {
    std::vector <int> seqSet;
    for(int i = 0; i < (int)seqList.list.size(); ++ i) {
        seqSet.emplace_back(i);
    }

    treeInBFS = std::queue <SequenceClusterNode> ();
    treeInBFS.push(SequenceClusterNode(seqList, seqSet, 0, 0.0));
}


SequenceClusterTree::~SequenceClusterTree() { }