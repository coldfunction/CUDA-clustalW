#ifndef PAIRWISEALIGNINFO
#define PAIRWISEALIGNINFO
#include <cstring>
#include <iostream>

using namespace clustalw;
using namespace std;


class pairwiseAlignInfoStruct{
public:
    int numSeqs;
    float pwGapOpen;
    float pwGapExtend;
    int gapPos1;
    int gapPos2;
    int matAvgScore;
    int DNAFlag;
    PairScaleValues scaleValues;


    pairwiseAlignInfoStruct(){
    }

    pairwiseAlignInfoStruct(
        int numseq,
        float pwgopen,
        float pwgextend,
        int gpos1,
        int gpos2,
        int mavgscore,
        PairScaleValues scale,
        int dnaflag):
              numSeqs(numseq), 
              pwGapOpen(pwgopen), 
              pwGapExtend(pwgextend), 
              gapPos1(gpos1), 
              gapPos2(gpos2),
              matAvgScore(mavgscore),
              scaleValues(scale),
              DNAFlag(dnaflag)
    {
    
    }

    void print_pair(){
        cout << "numSeqs:" <<  numSeqs << endl;
        cout << "pwGapOpen:" << pwGapOpen << endl;
        cout << "pwGapExtend:" << pwGapExtend << endl;
        cout << "gapPos1:" << gapPos1 << endl;
        cout << "gapPos2:" << gapPos2 << endl;
        cout << "matAvgScore:" << matAvgScore << endl;
        cout << "DNAFlag:" << DNAFlag << endl;
        cout << "intScale:" << scaleValues.intScale << endl;
        cout << "gapOpenScale:" << scaleValues.gapOpenScale << endl;
        cout << "gapExtendScale:" << scaleValues.gapExtendScale << endl;

    }
 
};


class pairwiseAlignInfo{
public:
    pairwiseAlignInfoStruct* alignInfo;
    int isGPU;


    //CPU constructor
    pairwiseAlignInfo():isGPU(0), alignInfo(NULL){
    
    }
 
    pairwiseAlignInfo(const pairwiseAlignInfoStruct& copy){
       isGPU = 0;
       alignInfo = new pairwiseAlignInfoStruct;
       memcpy(alignInfo, &copy, sizeof(pairwiseAlignInfoStruct));
    }
 
    pairwiseAlignInfo(const pairwiseAlignInfo& copy){
        isGPU = 0;
        alignInfo = new pairwiseAlignInfoStruct;
        memcpy(alignInfo, copy.alignInfo, sizeof(pairwiseAlignInfo));
    }

    //GPU constructor
    pairwiseAlignInfo(pairwiseAlignInfoStruct* CPUInfo, int GPUORNOT){
        if(GPUORNOT==1){
            isGPU = 1;
            cudaMalloc((void**)& alignInfo, sizeof(pairwiseAlignInfoStruct) );
            cudaMemcpy(alignInfo, CPUInfo, sizeof(pairwiseAlignInfoStruct), cudaMemcpyHostToDevice);
        }
        else    std::cout << "Can not allocate GPU memory. \n";
    }

    pairwiseAlignInfo(pairwiseAlignInfo* CPUInfo, int GPUORNOT) {
        if(GPUORNOT == 1){
            isGPU = 1;
            cudaMalloc((void**) &alignInfo, sizeof(pairwiseAlignInfoStruct));
            cudaMemcpy(alignInfo, CPUInfo->alignInfo, sizeof(pairwiseAlignInfoStruct), cudaMemcpyHostToDevice);
        }
        else    std::cout << "Can not allocate GPU memory. \n";
    }

    ~pairwiseAlignInfo(){
        if(isGPU == 1){
            cudaFree(alignInfo);
        }
        else{
            if(alignInfo!= NULL){
                delete alignInfo;
                alignInfo = NULL;
            }
        }
    }

};


#endif

