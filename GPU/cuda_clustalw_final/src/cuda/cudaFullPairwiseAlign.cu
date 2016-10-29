#include "../alignment/Alignment.h"
#include "../general/clustalw.h"
#include "cudaFullPairwiseAlign.h"
#include "DyArray2D.h"
#include "DyArray1D.h"
#include "pairwiseAlignInfo.h"
#include "Stack.cu"

#include <cuda.h>
#include <omp.h>
#include <iostream>
#include <vector>

#define SPThread 0
#define MAXLENGTH 1600
//#define MAXLENGTH 7950

#define numberOfBlocks  128 
#define numberOfThreads 32
#define nstreams        4


#define Matrix(k) tex1Dfetch(texMatrix, k)
texture<int, 1, cudaReadModeElementType> texMatrix;

__device__ __constant__ pairwiseAlignInfoStruct gpair;
using namespace std;


/* CPU function */
int findMaxLength(vector<vector <int> >& vec ){
    int maxLen = 0;
    for(int i = 0 ; i < vec.size(); i ++ ){
        if(vec[i].size() > maxLen)
            maxLen = vec[i].size();
    }
    return maxLen ; 
}

void Set_Element(int2* element, int element_size, int numSeqs){
    int index = 0;
    for(int i = 0 ; i<numSeqs ; i++){
        for(int j = i + 1 ; j<numSeqs ; j++){
            element[index] = make_int2(i, j);
            index++;
        }
    }
}



/* for single thread use */
__device__ int mmLength(unsigned char* Sequence,int length,  int gapPos1, int gapPos2){
    int len = 0;
    for (int i = 1; i <= length; i++){
        if ((Sequence[i] != gapPos1) && (Sequence[i] != gapPos2)) {
            len ++;
        }
    }
    return len;
}

__device__ int2 getGapOpenGapExtend(int DNAFlag, int matAvgScore, PairScaleValues scaleValues, float pwGapOpen, float pwGapExtend, int m, int n){
    int intScale       = scaleValues.intScale;
    int gapOpenScale   = scaleValues.gapOpenScale;
    int gapExtendScale = scaleValues.gapExtendScale;
    
    if (DNAFlag){
        return make_int2( static_cast<int>(2 * pwGapOpen * intScale * gapOpenScale),  static_cast<int>(pwGapExtend * intScale * gapExtendScale));
    }
    else {
        int gapOpen = 0;
        if (matAvgScore <= 0)
            gapOpen = 2 * static_cast<int>((pwGapOpen +__logf(static_cast<float>(min(n,m)))) * intScale);
        else
            gapOpen = static_cast<int>(2 * matAvgScore * (pwGapOpen + __logf(static_cast<float>(min(n,m)))) * gapOpenScale);

        return make_int2(gapOpen , static_cast<int>(pwGapExtend * intScale));
    }
}

__device__ int calcScore(unsigned char* sequence1, unsigned char* sequence2, int i, int j){
    return Matrix(sequence1[i] * 32 + sequence2[j]);
}

__device__ int cudaCalcScore(unsigned char* sequence1, unsigned char* sequence2, int iat, int jat, int v1, int v2){
    return calcScore(sequence1, sequence2, v1 + iat, v2 + jat);
}

__device__ void cudaAdd(int& lastPrint, int& printPtr, int* displ, int v){
    if (lastPrint < 0){
        displ[printPtr - 1] = v;
        displ[printPtr] = lastPrint;
    }
    else{
        lastPrint = displ[printPtr] = v;
    }
    printPtr++;
}

__device__ void cudaDel(int& lastPrint, int& printPtr, int* displ, int k){
    if (lastPrint < 0){
        lastPrint = displ[printPtr - 1] -= k;
    }
    else{
        lastPrint = displ[printPtr] =  - (k);
        printPtr++;
    }
}

__device__ int Gap(int k, int constant , int gapExtend){
    if(k <= 0)     return 0;
    else           return constant + gapExtend * k;
}
__device__ int cudaGap(int k, int gapOpen, int gapExtend){
    return Gap(k, gapOpen, gapExtend);
}

__device__ int cudaTbGap(int k, int tb, int gapExtend){
    return Gap(k, tb, gapExtend);
}

__device__ int cudaTeGap(int k, int te, int gapExtend){
    return Gap(k, te, gapExtend);
}

__device__ float cudaTracePath(unsigned char* sequence1, unsigned char* sequence2, int sb1, int sb2, int printPtr, int* displ, int gapPos1, int gapPos2){
    int toDo = printPtr - 1;
    int i1 = sb1;
    int i2 = sb2;

    int count = 0;
    for (int i = 1 ; i <= toDo ; i++){
        int k = displ[i];
        if (k == 0){
            int res1 = sequence1[i1];
            int res2 = sequence2[i2];

            if ((res1 != gapPos1) && (res2 != gapPos2) && (res1 == res2)){
                count++;
            }
            ++i1;
            ++i2;
        }
        else{
            if (k > 0)     i2 += k;
            else           i1 -= k;
        }
    }
    return 100.f*(float)count;
}


/* for multi thread function*/

__device__  void setSequence(unsigned char* Sequence, int* array2D, int length){
    int numThreads = blockDim.x;
    int threadId = threadIdx.x;

    for(int i = threadId + 1;  i <= length ; i += numThreads){
        Sequence[ i ] = static_cast<unsigned char> (array2D[ i]);
    }
}

__device__ void cudaInitHDOrRS(int* HR, int* DS, int gapOpen, int gapExtend, int tbte, int iStart, int iEnd, bool isInc){
    int numThreads = blockDim.x;
    int threadId = threadIdx.x;

    int k = (isInc ? 1 : -1);
    int i = iStart + threadId*k;
    
    if(threadId == SPThread )    HR[iStart-k] = 0;
    int t =  - tbte;
    while((isInc && i<=iEnd)|| (!isInc && i>=iEnd)){
        int value = (isInc ? t - i*gapExtend : t - (iStart - i + 1)*gapExtend );
        
        HR[i] = value;
        DS[i] = value - gapOpen;
        i += numThreads*k;
    }
}

__device__ void cudaSetHD(int* HH, int* DD, unsigned char* sequence1, unsigned char* sequence2, int gapOpen, int gapExtend, int tb, int A, int B, int iStart, int iEnd, int jStart, int jEnd){
    int numThreads = blockDim.x;
    int threadId = threadIdx.x;

    __shared__ bool over_loop;
    __shared__ bool go;


    if(threadId == SPThread){
        over_loop = false;
        go = false;
        if(iStart > iEnd)    over_loop = true;
    }

    int temp_count = ( threadId + iStart > iEnd ? 0 : threadId + 1 );
    int t =  - tb;

    __syncthreads();
    for (int i = threadId + iStart ; !over_loop ; i += numThreads) {
        int s = 0;
        int value = 0;
        int hh = 0;
        int  f = 0;
        int  g = 0;
        int  e = 0;
        
        int count = temp_count;
        bool writeHH0 = false;
        bool start = false;
        bool last = false;
        for (int j = jStart ; !go ; j++){
            count --;

            if(i <= iEnd && count == 0){
                start = true;
                if(i<iEnd && threadId == numThreads -1)     last = true;
                if(i == iEnd)                               last = true;
            }
            else if(count > 0 && i<= iEnd){
                start = false;
                j = jStart - 1;
            }
            
            if(start){
               if(!writeHH0){
                    s = HH[jStart - 1];
                    value = t - i*gapExtend;
                    hh = value;
                    f = value - gapOpen;
                    HH[jStart - 1] = hh;
                    writeHH0 = true;
               }
               
               f -= gapExtend;
               g = hh - gapOpen - gapExtend;
               if(f < g)    f = g;

               e = DD[j] - gapExtend;
               g = HH[j] - gapOpen - gapExtend;
               if(e < g)   e = g;
               
               hh = s + cudaCalcScore(sequence1, sequence2, i, j, A, B);
               
               if (f > hh)    hh = f;
               if (e > hh)    hh = e;

               s = HH[j];
               HH[j] = hh;
               DD[j] = e;

            }
            
            __syncthreads();
            if(j==jEnd   && i <= iEnd){
                start = false;
                if(last){
                    go = true;
                    if(i==iEnd)    over_loop = true;
                }
            }
            __syncthreads();
        }
        __syncthreads();
        if(threadId == SPThread){
            go = false;
        }
        __syncthreads();
    }
    if(threadId == SPThread)    DD[jStart - 1] = HH[jStart - 1];
    __syncthreads();

}


__device__ void cudaSetRS(int* RR, int* SS, unsigned char* sequence1, unsigned char* sequence2, int gapOpen, int gapExtend, int te, int A, int B, int iStart, int iEnd, int jStart, int jEnd){

    int numThreads = blockDim.x;
    int threadId = threadIdx.x;

    __shared__ bool over_loop;
    __shared__ bool go;

    if(threadId == SPThread){
        over_loop = false;
        go = false;
        if(iStart < iEnd)    over_loop = true;
    }

    int temp_count = (iStart - threadId >=  iEnd ? threadId + 1 : 0);
    int t =  - te;

    __syncthreads();
    for (int i = iStart - threadId ; !over_loop ; i -= numThreads) {
        int s =  0; 
        int value = 0;
        int hh = 0;
        int  f = 0;
        int  g = 0; 
        int  e = 0;

        int count = temp_count;
        bool writeRRN = false;
        bool start = false;
        bool last = false;
        for (int j = jStart ; !go ; j--){
            count --;

            if(i >= iEnd && count == 0){
                start = true;
                if( i>iEnd && threadId == numThreads -1)    last = true;
                if(i == iEnd)                               last = true;
            }
            else if(count > 0 && i>= iEnd){
                start = false;
                j = jStart + 1;
            }
            
            if(start){
                if(!writeRRN){
                    s = RR[jStart+1];
                    value = t - (iStart - i + 1)*gapExtend;
                    hh = value;
                    f = value - gapOpen;
                    RR[jStart + 1] = hh;
                    writeRRN = true;
                }
                
                f -= gapExtend;
                g = hh - gapOpen - gapExtend;
                if(f < g)   f = g;

                e = SS[j] - gapExtend;
                g = RR[j] - gapOpen - gapExtend;
                if (e < g)  e = g;
                
                hh = s + cudaCalcScore(sequence1, sequence2, i + 1, j + 1, A, B);
                if (f > hh)    hh = f;
                if (e > hh)    hh = e;

                s = RR[j];
                RR[j] = hh;
                SS[j] = e;
            }
            
            
            __syncthreads();
            if(j==jEnd){
                start = false;
                if(last){
                    go = true;
                    if(i==iEnd)    over_loop = true;
                }
            }
            __syncthreads();
        }
        __syncthreads();
        if(threadId == SPThread){
            go = false;
        }
        __syncthreads();
    }
    
    if(threadId == SPThread)    SS[jStart+1] = RR[jStart+1];
    __syncthreads();
}

__device__ int3 cudaForwardPath(unsigned char* sequence1, unsigned char* sequence2, int n, int m, int gapOpen, int gapExtend, int* HH, int* DD){
    int numThreads = blockDim.x;
    int threadId = threadIdx.x;

    int maxScore = 0;
    int se1 = 0;
    int se2 = 0;
    __shared__ bool over_loop;
    __shared__ bool go;


    __shared__ int3 S_ScoreAndPos;

    if(threadId == SPThread ){
        S_ScoreAndPos = make_int3(0, 0, 0);
        over_loop = false;
        go = false;
        if(n<1)    over_loop = true;
    }

    for(int i = threadId ;  i <= m ; i += numThreads){
        HH[i] = 0;
        DD[i] = -gapOpen;
    }
    

    int temp_count = ( threadId >= n ? 0 : threadId + 1 );
    
    __syncthreads();
    for (int i = threadId + 1 ; !over_loop ; i+= numThreads){
        int hh = 0;
        int  p = 0;
        int  f =  -gapOpen;
        int  t = 0;
        bool start = false;
        bool last = false;
        int count = temp_count;
        for (int j = 1  ; !go  ; j++) {
            count --;

            if(i <= n && count == 0){//can work
                start=true;
                if(i<n && threadId == numThreads -1)    last = true;
                if(i==n)                                last = true;
            }     
            else if(i<= n && count > 0){
                start = false;
                j = 0;
            }
            
            if(start){
                f -= gapExtend;
                t = hh - gapOpen - gapExtend;
                if (f < t)                       f = t;

                DD[j] -= gapExtend;
                t = HH[j] - gapOpen - gapExtend;
                if (DD[j] < t)                   DD[j] = t;

                hh = p + calcScore(sequence1, sequence2, i, j);
                
                if (hh < f)                      hh = f;
                if (hh < DD[j])                  hh = DD[j];
                if (hh < 0)                      hh = 0;

                p = HH[j];
                HH[j] = hh;

                if (hh > maxScore){
                    maxScore = hh;
                    se1 = i;
                    se2 = j;
                }
            }      
           
            __syncthreads();
            if(j==m && i <= n){
                start = false;
                if(maxScore > S_ScoreAndPos.x){
                    S_ScoreAndPos = make_int3(maxScore, se1, se2);
                }

                if(last){
                    go = true;
                    if(i == n)    over_loop = true;
                }
            }
            __syncthreads();
        }
        
        __syncthreads();
        if(threadId == SPThread){
            go = false;
        }
        __syncthreads();
    }
    
    __syncthreads();
    return S_ScoreAndPos;
}

__device__ int3 cudaReversePath(unsigned char* sequence1, unsigned char* sequence2, int n, int m, int gapOpen, int gapExtend, int* HH, int* DD, int3 ScoreAndPosSe){
    int numThreads = blockDim.x;
    int threadId = threadIdx.x;
    
    int maxScore = ScoreAndPosSe.x;
    int se1 = ScoreAndPosSe.y;
    int se2 = ScoreAndPosSe.z;

    int cost = 0;
    int sb1 = 1;
    int sb2 = 1;
    __shared__ bool over_loop;
    __shared__ bool go;
    __shared__ int3 S_ScoreAndPos;

    if(threadId == SPThread ){
        over_loop = false;
        go = false;
        S_ScoreAndPos = make_int3(0, 1, 1);
        if(se1 < 1)   over_loop = true;
    }

    for (int i = se2 - threadId ; i >= 1; i-= numThreads) {
        HH[i] = -1;
        DD[i] = -1;
    }

    int temp_count = (threadId >= se1 ? 0 : threadId+1);
    
    __syncthreads();
    for (int i = se1 - threadId ; !over_loop ; i-= numThreads) {

        int hh = -1;
        int  f = -1;
        int  p = (i == se1) ? 0 : -1;
        int  t = 0;

        bool start = false;
        bool last = false;
        int count = temp_count;
        
        for (int j = se2 ; !go ; j--) {
            count --;

            if(i>=1 && count == 0){
                start = true;
                if(i> 1 && threadId == numThreads -1)    last = true;
                if(i==1)                                 last = true;
            }
            else if(count > 0 && i >= 1){
                start = false;
                j = se2+1;
            }
            
            if(start){
                f -= gapExtend;
                t = hh - gapOpen - gapExtend;
                if (f < t)                           f = t;

                DD[j] -= gapExtend;
                t = HH[j] - gapOpen - gapExtend;
                if (DD[j] < t)                       DD[j] = t;

                hh = p + calcScore(sequence1, sequence2, i, j);
                if (hh < f)                          hh = f;
                if (hh < DD[j])                      hh = DD[j];

                p = HH[j];
                HH[j] = hh;

                if (hh > cost) {
                    cost = hh;
                    sb1 = i;
                    sb2 = j;
                    if (cost >= maxScore)    start = false;
                }
            }
            
            __syncthreads();
            if(j==1 && i>= 1){
                start = false;
                if(cost > S_ScoreAndPos.x){
                    S_ScoreAndPos = make_int3(cost, sb1, sb2);
                }
                
                if(last){
                    go = true;
                    if(i==1)    over_loop = true;
                }
            }
            __syncthreads();
        }
        
        __syncthreads();
        if(threadId == SPThread){
            go = false;
           // if (cost >= maxScore)    go = true;
        }
        __syncthreads();
    }
    __syncthreads();

    return S_ScoreAndPos;
}


__device__ 
int cudaDiff(unsigned char* sequence1, unsigned char* sequence2, int* HH, int* DD, int* RR, int* SS, int* displ, int gapOpen, int gapExtend, int3 ScoreAndPosSe, int3 ScoreAndPosSb){
    int threadId = threadIdx.x;
    __shared__ Stack<diffArgv> stack;
    __shared__ bool go;
    __shared__ int lastPrint;
    __shared__ int printPtr;
    __shared__ bool type;
    __shared__ diffArgv argument;
    __shared__ int3 ScoreAndPosMid;

    int& A = argument.A;
    int& B = argument.B;
    int& M = argument.M;
    int& N = argument.N;
    int& tb = argument.tb;
    int& te = argument.te;
    int& control = argument.isDel2;
 
    int midi = 0 ;
    int midj = 0;
    int midh = 0;
    int hh = 0;
    
    if(threadId == SPThread){
        diffArgv initArgv = {ScoreAndPosSb.y - 1, ScoreAndPosSb.z - 1, ScoreAndPosSe.y - ScoreAndPosSb.y + 1, ScoreAndPosSe.z - ScoreAndPosSb.z + 1, 0, 0, 0};
        StackInit(&stack);
        StackPush(&stack, initArgv);
        lastPrint = 0;
        printPtr = 1;
        go = false;
        type = true;
        ScoreAndPosMid = make_int3(0, 0, 0);
    }
    __syncthreads();
    while(stack.top != -1){
        __syncthreads();
        if(threadId == SPThread){
             argument = StackPop(&stack);
             type = true;
             go = false;
             
             if(control)      cudaDel(lastPrint, printPtr, displ, 2);
        
             if(N <= 0){
                 if(M > 0)    cudaDel(lastPrint, printPtr, displ, M);
                 go = true;
             }
             else if(M <= 1){
                 if(M <= 0)    cudaAdd(lastPrint, printPtr, displ, N);
                 else{
                    midh =  - (tb + gapExtend) - cudaTeGap(N, te, gapExtend);
                    hh   =  - (te + gapExtend) - cudaTbGap(N, tb, gapExtend);
                    if (hh > midh)     midh = hh; 
                    
                    midj = 0;
                    for (int j = 1; j <= N; j++) {
                        hh = cudaCalcScore(sequence1, sequence2, 1, j, A, B) 
                             - cudaTeGap(N - j, te, gapExtend) 
                             - cudaTbGap(j - 1, tb, gapExtend);

                        if (hh > midh) {
                            midh = hh;
                            midj = j;
                        }
                    }

                    if (midj == 0) {
                        cudaDel(lastPrint, printPtr, displ, 1);
                        cudaAdd(lastPrint, printPtr, displ, N);
                    }
                    else{
                        if(midj > 1)    cudaAdd(lastPrint, printPtr, displ, midj - 1);
                        
                        displ[printPtr++] = lastPrint = 0;
                        
                        if(midj < N)    cudaAdd(lastPrint, printPtr, displ, N - midj);
                    }
                 }
                 go = true;
             }
             
             ScoreAndPosMid.x = midh;
             ScoreAndPosMid.z = midj;
        }
               
        __syncthreads();
        midh = ScoreAndPosMid.x;
        midj = ScoreAndPosMid.z;
        
       
        if(go)     continue;

       
        //initialize
        cudaInitHDOrRS(HH, DD, gapOpen, gapExtend, tb, 1, N, true);
        cudaInitHDOrRS(RR, SS, gapOpen, gapExtend, te, N-1, 0, false);

        midi = M / 2;
       
        //Setting 
        cudaSetHD(HH, DD, sequence1, sequence2, gapOpen, gapExtend, tb, A, B,     1, midi,     1, N);
        cudaSetRS(RR, SS, sequence1, sequence2, gapOpen, gapExtend, te, A, B, M - 1, midi, N - 1, 0);
       
 
        //doing 
        midh = HH[0] + RR[0];
        midj = 0;
        if(threadId == SPThread ){
            for (int j = 0; j <= N; j++){
                hh = HH[j] + RR[j];
                if (hh >= midh){
                    if (hh > midh || (HH[j] != DD[j] && RR[j] == SS[j])){
                        midh = hh;
                        midj = j;
                    }
                }
            }

            for (int j = N ; j >= 0; j--){
                hh = DD[j] + SS[j] + gapOpen;
                if (hh > midh){
                    midh = hh;
                    midj = j;
                    type = false;
                }
            }
            ScoreAndPosMid.x = midh;
            ScoreAndPosMid.z = midj;
        }

        __syncthreads();
        midh = ScoreAndPosMid.x;
        midj = ScoreAndPosMid.z;

        if(threadId==SPThread){
            if (type) {             // Type 1 gaps
                diffArgv rArgv = {A + midi, B + midj, M - midi, N - midj, gapOpen, te, 0};
                diffArgv lArgv = {A, B, midi, midj, tb, gapOpen, 0};

                StackPush(&stack, rArgv);
                StackPush(&stack, lArgv);
            }
            else {
                diffArgv rArgv = {A + midi + 1, B + midj, M - midi - 1, N - midj, 0, te, 1};
                diffArgv lArgv = {A, B, midi - 1, midj, tb, 0, 0};

                StackPush(&stack, rArgv);
                StackPush(&stack, lArgv);
            }
        }
        __syncthreads();
    }
    __syncthreads();
    return printPtr;
}



__global__ void cudaFullPairwiseAlignKernel( 
    DyArray2DStruct<int>* sequence,
    float* distMatrix,
    DyArray2DStruct<int>* HH, 
    DyArray2DStruct<int>* DD, 
    DyArray2DStruct<int>* RR, 
    DyArray2DStruct<int>* SS, 
    DyArray2DStruct<int>* displ, 
    int2* element)
{
    int threadId = threadIdx.x;
    int blockId = blockIdx.x;
    
    __shared__ int2 gapOpen_gapExtend;
    int& gapOpen = gapOpen_gapExtend.x; // scaled to be an integer, this is not a mistake
    int& gapExtend = gapOpen_gapExtend.y; // scaled to be an integer, not a mistake
 
    int len1 = 0;
    int len2 = 0;
    int2 which = element[blockId];
    int si = which.x;
    int sj = which.y;
   
    int3 ScoreAndPosSe = make_int3(0, 0, 0);
    int3 ScoreAndPosSb = make_int3(0, 1, 1);

    int width = sequence->width;
    int indexHDRS =  HH->width * (blockId);
    int indexDispl = displ->width*(blockId);
 
    __shared__  bool control;

    __shared__ unsigned char sequence1[MAXLENGTH];
    __shared__ unsigned char sequence2[MAXLENGTH];
    

    int n = sequence->array2DSize[si + 1] - 1;
    int index = (si + 1)*width;
    setSequence(sequence1, &sequence->array2D[index],  n);

    // copy sequence2
    int m = sequence->array2DSize[sj + 1] - 1;
    index = (sj + 1)*width;
    setSequence(sequence2, &sequence->array2D[index], m);

    __syncthreads();
    if(threadId == SPThread){
        len1 = mmLength(sequence1, n, gpair.gapPos1, gpair.gapPos2);
        len2 = mmLength(sequence2, m, gpair.gapPos1, gpair.gapPos2);
        control = false;
        if (n == 0 || m == 0){
            distMatrix[blockId] = 1.0f;
            control = true;
        }
        gapOpen_gapExtend = getGapOpenGapExtend(gpair.DNAFlag, gpair.matAvgScore, gpair.scaleValues, gpair.pwGapOpen, gpair.pwGapExtend, m, n);
    }
    __syncthreads();
            
    if(control)    return;

    ScoreAndPosSe = cudaForwardPath(sequence1, sequence2, n, m, gapOpen, gapExtend, &HH->array2D[indexHDRS],& DD->array2D[indexHDRS]); 

    ScoreAndPosSb = cudaReversePath(sequence1, sequence2, n, m, gapOpen, gapExtend, &HH->array2D[indexHDRS],& DD->array2D[indexHDRS], ScoreAndPosSe); 
    int printPtr = cudaDiff(
                       sequence1,
                       sequence2, 
                       &HH->array2D[indexHDRS], 
                       &DD->array2D[indexHDRS], 
                       &RR->array2D[indexHDRS], 
                       &SS->array2D[indexHDRS], 
                       &displ->array2D[indexDispl],
                       gapOpen, gapExtend, 
                       ScoreAndPosSe, ScoreAndPosSb);

    if(threadId == SPThread){
        double mmScore = (double)cudaTracePath(sequence1, sequence2, ScoreAndPosSb.y, ScoreAndPosSb.z, printPtr, &displ->array2D[indexDispl], gpair.gapPos1, gpair.gapPos2);
        if ((len1 == 0) || (len2 == 0))    mmScore = 0.f;
        else                               mmScore = 100 - __fdiv_rn (mmScore, static_cast<float>(min(len1, len2)));
        
        float score = mmScore/100.f;
        distMatrix[blockId] = score;
    }

}


void cudaFullPairwiseAlign(Alignment* alignPtr, DistMatrix* distMat, int iStart , int iEnd, int jStart, int jEnd){
    if(distMat->getSize() != alignPtr->getNumSeqs() + 1){
        printf( "The distance matrix is not the right size!\n");
        printf("Need to terminate program.\n");
        exit(1);
    }
    if((iStart < 0) || (iEnd < iStart) || (jStart < 0) || (jEnd < jStart)){
        printf("The range for pairwise Alignment is incorrect.\n");
        printf("Need to terminate program.\n");
        exit(1);
    }
   
    if(alignPtr->getNumSeqs() == 0){
        return;
    }

    int matrix[NUMRES][NUMRES];
    int _matAvgScore;
    PairScaleValues scaleValues;


    int maxRes = subMatrix->getPairwiseMatrix(matrix, scaleValues, _matAvgScore);
    if (maxRes == 0){
        printf("Could not get the substitution matrix\n");
        return;
    }
 
    int maxAlnLength = alignPtr -> getMaxAlnLength(); 
    int numSeqs = alignPtr->getNumSeqs();

    int element_size = numSeqs*(numSeqs-1)/2;
    int2* element = new int2 [element_size];
    Set_Element(element, element_size, numSeqs);


    int create_num_threads = 0;
    int num_device = 0;
    cudaGetDeviceCount(&num_device);

    if(num_device < 1){
        printf("no CUDA capable devices were detected\n");
        exit(0);
    }

    cout << "Run " << num_device << " device" << endl;
    create_num_threads = num_device;
    omp_set_num_threads(create_num_threads);
    #pragma omp parallel
    {

    int threadid = omp_get_thread_num();
   // cout << "threadid: "<< threadid << endl;

    cudaError_t error2 = cudaSetDevice(threadid % num_device);
   // cout << "set device: " << cudaGetErrorString(error2) << endl;
    int gpu_id = -1;
    cudaGetDevice(&gpu_id);
    if(gpu_id != (threadid %num_device)){
        printf("Set device error \n");
        exit(0);
    }
    //streams
    cudaStream_t stream[nstreams];
    for(int k = 0 ; k < nstreams ; k++){
        if(cudaStreamCreate(&stream[k]) != cudaSuccess){
            cout << "Error \n";
        }
    }

    
    //create openmp thread
    //texture for matrix[32][32]
    int* gMatrix;
    cudaMalloc((void**)&gMatrix, sizeof(int)*NUMRES*NUMRES );
    for(int i = 0 ; i < 32 ; i += 32)
        cudaMemcpy(&gMatrix[i], &matrix[i], sizeof(int)*NUMRES*NUMRES, cudaMemcpyHostToDevice);
    cudaBindTexture(0, texMatrix, gMatrix, NUMRES*NUMRES*sizeof(int));

   
   
    //pair  constant memory
    pairwiseAlignInfoStruct pair(
            alignPtr->getNumSeqs(),
            userParameters->getPWGapOpen(),
            userParameters->getPWGapExtend(),
            userParameters->getGapPos1(),  
            userParameters->getGapPos2(),
            _matAvgScore,
            scaleValues,
            (int)(userParameters->getDNAFlag()));
    cudaMemcpyToSymbol((char*)&gpair, &pair, sizeof(pairwiseAlignInfoStruct));


    //sequence
    vector<vector<int> >& ptrToSequence = *const_cast< vector<vector<int> >* >(alignPtr->getSeqArray());
    DyArray2DStruct<int> sequenceConvertToArray( &ptrToSequence, findMaxLength(ptrToSequence));
    DyArray2DStruct<int> sequence(sequenceConvertToArray , 1); 
    DyArray2D<int> gSequence(&sequence, 1);



    //distMatrix
    //create pin-lock memory
    float* distMatrix;
    cudaMallocHost((void**)&distMatrix, sizeof(float)*numberOfBlocks*nstreams); 
    memset(distMatrix, 0, sizeof(float)*numberOfBlocks*nstreams);

    float* gDistMatrix;
    cudaError_t error = cudaMalloc((void**)&gDistMatrix, sizeof(float)*numberOfBlocks*nstreams);
    cudaMemset(gDistMatrix, 0, sizeof(float)*numberOfBlocks*nstreams);
    
    // pairwise needs global memory in every blocks
    DyArray2DStruct<int>* HHarray = new DyArray2DStruct<int>[nstreams](); 
    DyArray2DStruct<int>* DDarray = new DyArray2DStruct<int>[nstreams](); 
    DyArray2DStruct<int>* RRarray = new DyArray2DStruct<int>[nstreams]();
    DyArray2DStruct<int>* SSarray = new DyArray2DStruct<int>[nstreams]();
    DyArray2DStruct<int>* displarray = new DyArray2DStruct<int>[nstreams]();

    
    for(int k = 0 ; k < nstreams ; k++){
        HHarray[k].initgpu (maxAlnLength, numberOfBlocks, 1);
        DDarray[k].initgpu (maxAlnLength, numberOfBlocks, 1); 
        RRarray[k].initgpu (maxAlnLength, numberOfBlocks, 1);
        SSarray[k].initgpu (maxAlnLength, numberOfBlocks, 1);
        displarray[k].initgpu ((2 * maxAlnLength) + 1, numberOfBlocks, 1);
    }

    DyArray2D<int>* gHHarray = new DyArray2D<int>[nstreams]();
    DyArray2D<int>* gDDarray = new DyArray2D<int>[nstreams]();
    DyArray2D<int>* gRRarray = new DyArray2D<int>[nstreams]();
    DyArray2D<int>* gSSarray = new DyArray2D<int>[nstreams]();
    DyArray2D<int>* gdisplarray = new DyArray2D<int>[nstreams]();

    for(int k = 0 ; k < nstreams ; k++){
        gHHarray[k].initgpu(&HHarray[k], 1);
        gDDarray[k].initgpu(&DDarray[k], 1);
        gRRarray[k].initgpu(&RRarray[k], 1);
        gSSarray[k].initgpu(&SSarray[k], 1);
        gdisplarray[k].initgpu(&displarray[k], 1);
    }

    int2* gElement;
    error = cudaMalloc((void**)&gElement, sizeof(int2)*element_size );
    error = cudaMemcpy(gElement, element, sizeof(int2)*element_size, cudaMemcpyHostToDevice);


    int numberOfStreamBlocks = nstreams * numberOfBlocks;

    cudaEvent_t start_event, stop_event;

    cudaEventCreate(&start_event);
    cudaEventCreate(&stop_event);

    cudaEventRecord(start_event, 0);


    #pragma omp for //schedule(dynamic, 1)
    for(int i = 0; i<element_size ; i+= numberOfStreamBlocks){
        int rest = element_size - i;
        int rest2 = element_size - i;

 
        for(int k = 0 ; k < nstreams ; k++){
            if(rest >= numberOfBlocks){
                cudaFullPairwiseAlignKernel<<<numberOfBlocks,  numberOfThreads, 0, stream[k]>>>(
                        gSequence.array2DPtr, 
                        &gDistMatrix[k*numberOfBlocks], 
                        gHHarray[k].array2DPtr, 
                        gDDarray[k].array2DPtr,
                        gRRarray[k].array2DPtr,
                        gSSarray[k].array2DPtr,
                        gdisplarray[k].array2DPtr, 
                        &gElement[i + k*numberOfBlocks]);
            }
            else{
                cudaFullPairwiseAlignKernel<<<rest,  numberOfThreads, 0, stream[k]>>>(
                        gSequence.array2DPtr, 
                        &gDistMatrix[k*numberOfBlocks], 
                        gHHarray[k].array2DPtr, 
                        gDDarray[k].array2DPtr,
                        gRRarray[k].array2DPtr,
                        gSSarray[k].array2DPtr,
                        gdisplarray[k].array2DPtr, 
                        &gElement[i + k*numberOfBlocks]);
                break;
            }
            rest -= numberOfBlocks;
        }

        for(int k = 0 ; k < nstreams ; k++){
            if(rest2 >= numberOfBlocks){
                cudaMemcpyAsync(&distMatrix[k*numberOfBlocks], &gDistMatrix[k*numberOfBlocks], sizeof(float)*numberOfBlocks, cudaMemcpyDeviceToHost, stream[k]);
                //cudaMemcpy(&distMatrix[k*numberOfBlocks], &gDistMatrix[k*numberOfBlocks], sizeof(float)*numberOfBlocks, cudaMemcpyDeviceToHost);
            }
            else{
                cudaMemcpyAsync(&distMatrix[k*numberOfBlocks], &gDistMatrix[k*numberOfBlocks], sizeof(float)*rest2, cudaMemcpyDeviceToHost, stream[k]);
                //cudaMemcpy(&distMatrix[k*numberOfBlocks], &gDistMatrix[k*numberOfBlocks], sizeof(float)*rest2, cudaMemcpyDeviceToHost);
                break;
            }
            rest2 -= numberOfBlocks;
        }
        
        for(int k = 0; k<nstreams ; k++){
            cudaStreamSynchronize(stream[k]);
            for(int j=i+k*numberOfBlocks ; j<(i+(k+1)*numberOfBlocks) && j<element_size ; j++){
                distMat->SetAt(element[j].x + 1, element[j].y + 1, distMatrix[j-i]);
                distMat->SetAt(element[j].y + 1, element[j].x + 1, distMatrix[j-i]);
                if(userParameters->getDisplayInfo())
                {
                    utilityObject->info("Sequences (%d:%d) Aligned. Score:  %d",
                            element[j].x + 1, element[j].y+1, (int)(100.f - distMatrix[j-i]*100.f));
                }
            }
        }

    }


    cudaThreadSynchronize();
    cudaUnbindTexture(texMatrix);
    cudaFree(gElement);
    cudaFree(gDistMatrix);
    delete [] HHarray;
    delete [] DDarray;
    delete [] RRarray;
    delete [] SSarray;
    delete [] displarray;

    delete [] gHHarray;
    delete [] gDDarray;
    delete [] gRRarray;
    delete [] gSSarray;
    delete [] gdisplarray;


    for(int k = 0 ; k < nstreams ; k++){
        cudaStreamDestroy(stream[k]);
    }

    cudaEventRecord(stop_event, 0);
    cudaEventSynchronize(stop_event);

    float time_memcpy = 0.f;
    cudaEventElapsedTime(&time_memcpy, start_event, stop_event);
   // printf("device %d  pairwise alignment execute :\t%.2f ms\n",threadid,  time_memcpy);


    cout << cudaGetErrorString(cudaGetLastError()) << endl;

    }//end of parallel
 
    delete [] element;

    if(userParameters->getDisplayInfo())
    {
        distMat->printArray();
    }
}


#undef Matrix
