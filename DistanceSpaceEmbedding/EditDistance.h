#ifndef __EditDistance_DistanceSpaceEmbedding_H__
#define __EditDistance_DistanceSpaceEmbedding_H__

#include "Sequence.h"

class EditDistance
{
    private:
        const SequenceList seqList;
    public:
        EditDistance(const SequenceList&);
        ~EditDistance();

        int getEditDistance(const std::string&, const std::string&);
        void calculate(SequencePair&);
};


#endif