#ifndef CUDAMULTISEQALIGN_H
#define CUDAMULTISEQALIGN_H

#include "../alignment/Alignment.h"
#include "../tree/AlignmentSteps.h"
#include "../general/clustalw.h"

using namespace clustalw;
using namespace std;
#define SPThread 0

int cudaMultiSeqAlign(Alignment* alignPtr, DistMatrix* distMat, vector<int>* seqWeight, AlignmentSteps* progSteps, int iStart);

#endif

