#ifndef SPARSEMATRIXMAKER_CLUSTERMON_H_
#define SPARSEMATRIXMAKER_CLUSTERMON_H_

#include <utility>
#include <vector>
#include <algorithm>
#include "SeqCode.h"

class HomoPair {
	public:
		std::unordered_map <std::string, double> seqPair;

		HomoPair();
		~HomoPair();

		static std::string indexPairToString(const ID_TYPE&, const ID_TYPE&);
		static std::pair <ID_TYPE, ID_TYPE> stringToIndexPair(const std::string&);
		void create(const ID_TYPE&, const ID_TYPE&, const double&);
};


class SparseMatrixMaker {
    protected:
        const SeqCode& table;

    public:
        std::vector<ID_TYPE> cluster;
        HomoPair result;

        SparseMatrixMaker(const SeqCode&);
        ~SparseMatrixMaker();

        virtual void convert(const std::string&);
        void outputAsAllPairs(const std::string&);
};


class LinclustConv : public SparseMatrixMaker {
    public:
        LinclustConv(const SeqCode&);
        ~LinclustConv();

        void convert(const std::string&);
};


class ALFATConv : public SparseMatrixMaker {
    public:
        ALFATConv(const SeqCode&);
        ~ALFATConv();

        void convert(const std::string&);
};


class LSHConv : public SparseMatrixMaker {
    public:
        LSHConv(const SeqCode&);
        ~LSHConv();

        void convert(const std::string&);
};


class OtherConv : public SparseMatrixMaker {
    //T.B.A.
};

#endif
