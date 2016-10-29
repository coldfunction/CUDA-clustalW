/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef NJTREE_H
#define NJTREE_H

#include <vector>
#include <iomanip>
#include "ClusterTreeAlgorithm.h"
#include "../general/userparams.h"

namespace clustalw
{

class NJTree : public ClusterTreeAlgorithm
{
    public:
        NJTree(): verbose(false){};
        virtual void generateTree(clustalw::PhyloTree* phyTree, clustalw::DistMatrix* distMat, clustalw::SeqInfo* seqInfo,
                                  ofstream* tree = 0);
        virtual void setVerbose(bool choice){verbose = choice;};
    private:
        vector<double> av;
        vector<int> tkill;
        bool verbose;
};

}
#endif
