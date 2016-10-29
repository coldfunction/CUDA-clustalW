#ifndef CUDAFULLPAIRWISWALIGN_H
#define CUDAFULLPAIRWISWALIGN_H
#include "../alignment/Alignment.h"
#include "../general/clustalw.h"

using namespace clustalw;
using namespace std;

class diffArgv{
public:
    int A;
    int B;
    int M;
    int N;
    int tb;
    int te;
    int isDel2;
};

void cudaFullPairwiseAlign(Alignment* alignPtr, DistMatrix* distMat, int iStart , int iEnd, int jStart, int jEnd);

#endif

