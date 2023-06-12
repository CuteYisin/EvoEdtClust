#include "Sequence.h"


Sequence::Sequence(const std::string& header, const std::string& sequence, const int& index):
    header(header),
    sequence(sequence),
    index(index) {
}


Sequence::~Sequence() {}


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
