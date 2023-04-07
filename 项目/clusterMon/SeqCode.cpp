#include "SeqCode.h"

SeqCode::SeqCode() {
    seqID = std::unordered_map <std::string, ID_TYPE> ();
}


SeqCode::~SeqCode() { }


void SeqCode::code(const std::string& LookupTableFile) {
    std::cout << "+++ Load LookupTable file from " << LookupTableFile << std::endl;
	std::cout << "*** Code sequence to ID ..." << std::endl;

    std::ifstream ifs(LookupTableFile);
    std::string header, ID;
    int progress = 0;
    if(ifs.is_open()) {
        while(!ifs.eof()) {
            std::getline(ifs, header, '\t');
            std::getline(ifs, ID);
            if(header[0] == '>') {
                seqID[header] = std::stoi(ID);

                progress ++;
                if(progress % 1000000 == 0) {
					std::cout << "*** " << progress << " sequences coded ..." << std::endl;
				}
            }
        }
        std::cout << "*** Finally, " << progress << " sequences coded." << std::endl;
    } else {
        std::cerr << "!!! Error: failed to open " << LookupTableFile << std::endl;
    }
}
