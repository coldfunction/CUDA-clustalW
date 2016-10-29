#ifndef MSAINFO_H
#define MSAINFO_H
#include <cstring>
#include <iostream>

using namespace clustalw;
using namespace std;


class MSAInfoStruct{
public:
    int numSeqs;
   // int matrix[NUMRES][NUMRES];
   // float pwGapOpen;
  //  float pwGapExtend;
    int gapPos1;
    int gapPos2;
 //   int matAvgScore;
    int DNAFlag;
 //   PairScaleValues scaleValues;
    int extraEndElemNum;
    int ENDALN;
    int gapOpen;
    int gapExtend;
    int negMatrixFlag;
    int structPenalties1;
    int structPenalties2;
    int useSS1Flag;
    int useSS2Flag;
    string :q


    MSAInfoStruct(){
    }

    MSAInfoStruct(
        int numseq,
        int mtx[NUMRES][NUMRES], 
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
        for(int i = 0 ; i< NUMRES ; i++){
            for(int j = 0 ; j<NUMRES ; j++){
                matrix[i][j] = mtx[i][j];
            }
        }
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

        cout << "matrix[][]" << endl;
        for(int i = 0 ; i< NUMRES ; i++){
            for(int j = 0 ; j<NUMRES ; j++){
               cout << matrix[i][j] << "\t";
            }
            cout << endl;
        }
    }


 
};


class MSAInfo{
public:
    MSAInfoStruct* alignInfo;
    int isGPU;


    //CPU constructor
    MSAInfo():isGPU(0), alignInfo(NULL){
    
    }
 
    MSAInfo(const MSAInfoStruct& copy){
       isGPU = 0;
       alignInfo = new MSAInfoStruct;
       memcpy(alignInfo, &copy, sizeof(MSAInfoStruct));
    }
 
    MSAInfo(const MSAInfo& copy){
        isGPU = 0;
        alignInfo = new MSAInfoStruct;
        memcpy(alignInfo, copy.alignInfo, sizeof(MSAInfo));
    }

    //GPU constructor
    MSAInfo(MSAInfoStruct* CPUInfo, int GPUORNOT){
        if(GPUORNOT==1){
            isGPU = 1;
            cudaMalloc((void**)& alignInfo, sizeof(MSAInfoStruct) );
            cudaMemcpy(alignInfo, CPUInfo, sizeof(MSAInfoStruct), cudaMemcpyHostToDevice);
        }
        else    std::cout << "Can not allocate GPU memory. \n";
    }

    MSAInfo(MSAInfo* CPUInfo, int GPUORNOT) {
        if(GPUORNOT == 1){
            isGPU = 1;
            cudaMalloc((void**) &alignInfo, sizeof(MSAInfoStruct));
            cudaMemcpyToSymbol(alignInfo, CPUInfo->alignInfo, sizeof(MSAInfoStruct), cudaMemcpyHostToDevice);
            //cudaMemcpy(alignInfo, CPUInfo->alignInfo, sizeof(MSAInfoStruct), cudaMemcpyHostToDevice);
        }
        else    std::cout << "Can not allocate GPU memory. \n";
    }

    ~MSAInfo(){
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

