#ifndef HOMOSEARCH_GROUNDTRUTH_H_
#define HOMOSEARCH_GROUNDTRUTH_H_

#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <vector>

typedef unsigned int ID_TYPE;

class HomoSearch {
    public:
        std::vector <std::string> seq;
        std::vector <std::pair <ID_TYPE, ID_TYPE>> HomoPair;

        HomoSearch();
        ~HomoSearch();

        void seqReservation(const std::string&);
        double EditSimilarity(const std::string&, const std::string&);
        void filter(const double&);
        void printResult(const std::string&);
};

#endif