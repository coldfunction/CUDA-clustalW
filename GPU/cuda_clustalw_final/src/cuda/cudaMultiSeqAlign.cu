#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#include "../general/clustalw.h"
#include "../alignment/Alignment.h"
#include "../tree/AlignmentSteps.h"
#include "../multipleAlign/MyersMillerProfileAlign.h"

#include "cudaMultiSeqAlign.h"
//#include "DyArray2D.h"
//#include "DyArray1D.h"
//#include "Stack.cu"

#include <cuda.h>
#include <iostream>
#include <vector>
#include <omp.h>
#include <cstdlib>

//#define numberOfThreads 1
//#define MAXLENGTH 1024
//#define numberOfBlocks 1

using namespace std;

void Random(vector<int>* v){
    vector<int>& V = *v;
    int size = V.size();

    srand(time(NULL));
    for(int i = 0 ; i < size ; i++){
        int idxa = rand()%size;
        int idxb = rand()%size;
        int temp = V[idxa];
        V[idxa]  = V[idxb];
        V[idxb] = temp;
    }
}

vector<int> Sort_Sets(vector< vector<int> > * ptrToSets){

    vector<vector <int> >& Sets = *ptrToSets;

    vector<int>    groupTable;
    vector<vector <int> > indexTable;
    int next = 1;
    int count = Sets.size();
    int size = count;
    int* flags = new int[size];
    memset(flags, 0, sizeof(int)*size);

    while(count > 1){
        int tracePath[Sets[1].size()];
        memset(tracePath, 0, sizeof(int)*Sets[1].size());
        for(int j  = 0 ; j < Sets[1].size() ; j ++){
            if(Sets[next][j] != 0)    tracePath[j] = 1;
        }

        vector<int> index;
        flags[next] = 1;
        index.push_back(next);
        count --;
        

        for(int k = 0 ; k < index.size() ; k ++){
            for(int i = 1 ; i< Sets.size() ; i ++){
                if(index[k] == i)    continue;
                if(flags[i] == 1)    continue;

                bool control = true;
                for(int j = 1 ; j < Sets[i].size() ; j++){
                    if((Sets[i][j]) != 0 && (tracePath[j]) != 0){
                        control = false;
                        break;
                    }
                }
                if(control){
                    for(int j = 1 ; j < Sets[i].size() ; j ++){
                        if(Sets[i][j] != 0){
                             tracePath [j] = 1;
                        }
                    }
                    index.push_back(i);
                    flags[i] = 1;
                    count --;
                }
            }
        }
        //Random(&index);
        indexTable.push_back(index);
        groupTable.push_back(index.size());

        for(int ii = 1 ; ii <  size; ii++){
            if(flags[ii] == 0){
                next = ii;
                break;
            }
        }
    }

    vector <vector<int> > buffer;
    vector<int> temp;
    buffer.push_back(temp);

    for(int i=0 ; i<indexTable.size() ; i++){
        for(int j = 0 ; j < indexTable[i].size() ; j ++){
            vector<int> temp = Sets[ indexTable[i][j] ];
            buffer.push_back(temp);
        }
    }

    Sets = buffer;

    delete [] flags;
    return groupTable;
}

/*
vector<int> Sort_Sets(vector< vector<int> > * ptrToSets){

    vector<vector <int> >& Sets = *ptrToSets;

    vector<int>    groupTable;
    vector<vector <int> > indexTable;
    int next = 1;
    int count = Sets.size();
    int size = count;
    int* flags = new int[size];
    memset(flags, 0, sizeof(int)*size);

    while(count > 1){
        vector<int> index;
        flags[next] = 1;
        index.push_back(next);
        count --;

        for(int k = 0 ; k < index.size() ; k ++){
            for(int i = 1 ; i< Sets.size() ; i ++){
                if(index[k] == i)    continue;
                if(flags[i] == 1)    continue;

                bool control = true;
                for(int j = 1 ; j < Sets[i].size() ; j++){
                    if((Sets[i][j]) != 0 && (Sets[ index[k] ][j]) != 0){
                        control = false;
                        break;
                    }
                }
                if(control){
                    // check inner
                    for(int inner = 0 ; inner < index.size() ; inner ++){
                        for(int j = 1 ; j < Sets[i].size() ; j ++){
                            if(Sets[i][j] != 0 && Sets[ index[inner] ][j] != 0){
                                control = false;
                                break;
                            }
                        }
                    }
                    if(control){
                        index.push_back(i);
                        flags[i] = 1;
                        count --;
                    }
                }
            }
        }
        indexTable.push_back(index);
        groupTable.push_back(index.size());

        for(int ii = 1 ; ii <  size; ii++){
            if(flags[ii] == 0){
                next = ii;
                break;
            }
        }
    }

    vector <vector<int> > buffer;
    vector<int> temp;
    buffer.push_back(temp);

    for(int i=0 ; i<indexTable.size() ; i++){
        for(int j = 0 ; j < indexTable[i].size() ; j ++){
            vector<int> temp = Sets[ indexTable[i][j] ];
            buffer.push_back(temp);
        }
    }

    Sets = buffer;

    delete [] flags;
    return groupTable;
}
*/

int cudaMultiSeqAlign(Alignment* alnPtr, DistMatrix* distMat, vector<int>* seqWeight, AlignmentSteps* progSteps, int iStart){

    if(!progSteps){
        return 0;
    }
    utilityObject->info("Start of Multiple Alignment\n"); 
        
    int* aligned;
    int ix;
    int numSeqs = alnPtr->getNumSeqs();
    int numSteps = progSteps->getNumSteps();
    
    int* maxid = new int[numSeqs + 1];
    vector<int> newOutputIndex(numSeqs);   
   
    alnPtr->addSeqWeight(seqWeight);
   // distMat->printArray(); 
    // for each sequence, find the most closely related sequence

    for (int i = 1; i <= numSeqs; i++){
        maxid[i] =  -1;
        for (int j = 1; j <= numSeqs; j++){
            if (j != i && maxid[i] < (*distMat)(i, j)){
                maxid[i] = static_cast<int>((*distMat)(i, j));
            }
        }
    }

    // group the sequences according to their relative divergence

    if (iStart == 0){

        // start the multiple alignments.........
        utilityObject->info("Aligning...");
        // first pass, align closely related sequences first....
        
        ix = 0;
        aligned = new int[numSeqs + 1];
        
        memset(aligned, 0, sizeof(int)*(numSeqs + 1));

        vector<vector<int> >* ptrToSets = const_cast<vector<vector<int> >* >(progSteps->getSteps());
        for (int set = 1 ; set <= numSteps; ++set) {
            for (int i = 1; i <= numSeqs; i++) {
                if (((*ptrToSets)[set][i] != 0) && (maxid[i] > userParameters->getDivergenceCutoff())) {
                    if (aligned[i] == 0) {
                        if (userParameters->getOutputOrder() == INPUT) {
                            ++ix;
                            newOutputIndex[i - 1] = i;
                        }
                        else {
                            if(ix >= newOutputIndex.size()) {
                                cout << "ERROR: size = " << newOutputIndex.size() << "ix = " << ix << "\n";
                                exit(1);
                            }
                            else {
                                newOutputIndex[ix] = i;
                                ++ix;
                            }
                        }
                        aligned[i] = 1;
                    }
                }
                if(aligned[i] == 0)    (*ptrToSets)[set][i] = 0; 
            }
        }

        
      //  progSteps->printAlignSteps();
        double sort_set_time = clock();
//        vector<int> groupTable = Sort_Sets(ptrToSets);
        sort_set_time = clock() - sort_set_time;
        cout << "sort_set_time: " << sort_set_time << endl;
       // progSteps->printAlignSteps();

        ProfileAlignAlgorithm* alignAlgorithm = new MyersMillerProfileAlign();
        int numProc = 1;
        //int numProc = omp_get_num_procs();
        //#pragma omp parallel num_threads(numProc) 
       // {
           // for ( int groupidx = 0, upbound = 1, start = 1 ; groupidx < groupTable.size() ; groupidx++) {
             //   upbound += groupTable[groupidx]; 
              
  //              #pragma omp for schedule(dynamic, 1)
               // for( int set = start ; set < upbound  ;  set ++){
                for(int set = 1 ; set <= numSteps ; set ++){
                    int entries = 0;
                    for (int i = 1; i <= numSeqs; i++) {
                        if (((*ptrToSets)[set][i] != 0) && (maxid[i] > userParameters->getDivergenceCutoff())) {
                            entries++;
                        }
                    }

                    int score = 0;
                    if (entries > 0) {
                        #if DEBUGFULL
                            if(logObject && DEBUGLOG) {    
                                logObject->logMsg("Doing profile align");
                            }
                        #endif
                        bool value = false;
//                        ProfileAlignAlgorithm* alignAlgorithm = new MyersMillerProfileAlign();
                        score = alignAlgorithm->profileAlign(alnPtr, distMat, progSteps->getStep(set), aligned);
//                        delete alignAlgorithm; 
                    }
                    else {
                        score = 0;
                    }

                    if(userParameters->getDisplayInfo()) {
                        if ((entries > 0) && (score > 0)) {
                            utilityObject->info("Group %d: Sequences:%4d      Score:%d", set, entries, score);                     
                        }
                        else {
                            utilityObject->info("Group %d:                     Delayed", set);
                        }
                    }
                }
       //         cout << endl;
         //       start = upbound;
           // }
       // }
          delete alignAlgorithm; 
    }
    else {
        aligned = new int[numSeqs + 1];
        ix = 0;
        for (int i = 1; i <= iStart + 1; i++) {
            aligned[i] = 1;
            ++ix;
            newOutputIndex[i - 1] = i;
        }
        for (int i = iStart + 2; i <= numSeqs; i++) {
            aligned[i] = 0;
        }
    }

    // second pass - align remaining, more divergent sequences..... 

    // if not all sequences were aligned, for each unaligned sequence,
    // find it's closest pair amongst the aligned sequences.

    vector<int> group;
    group.resize(numSeqs + 1); 
    vector<int> treeWeight;
    treeWeight.resize(numSeqs); 

    for (int i = 0; i < numSeqs; i++) {
        treeWeight[i] = (*seqWeight)[i];
    }

    // if we haven't aligned any sequences, in the first pass - align the
    // two most closely related sequences now
    if (ix == 0) {
        int max =  -1;
        int iseq = 0;

        for (int i = 1; i <= numSeqs; i++) {
            for (int j = i + 1; j <= numSeqs; j++) {
                if (max < (*distMat)(i, j)) {
                    max = static_cast<int>((*distMat)(i, j)); // Mark change 17-5-07
                    iseq = i;
                }
            }
        }
        aligned[iseq] = 1;
        if (userParameters->getOutputOrder() == INPUT) {
            ++ix;
            newOutputIndex[iseq - 1] = iseq;
        }
        else {
            newOutputIndex[ix] = iseq;
            ++ix;
        }
    }

    while (ix < numSeqs) {
        for (int i = 1; i <= numSeqs; i++) {
            if (aligned[i] == 0) {
                maxid[i] =  - 1;
                for (int j = 1; j <= numSeqs; j++){
                    if ((maxid[i] < (*distMat)(i, j)) && (aligned[j] != 0)) {
                        maxid[i] = static_cast<int>((*distMat)(i, j));// Mark change 17-5-07
                    }
                }
            }
        }
        // find the most closely related sequence to those already aligned

        int max =  - 1;
        int iseq = 0;
        for (int i = 1; i <= numSeqs; i++) {
            if ((aligned[i] == 0) && (maxid[i] > max)) {
                max = maxid[i];
                iseq = i;
            }
        }


        // align this sequence to the existing alignment
        // weight sequences with percent identity with profile
        // OR...., multiply sequence weights from tree by percent identity with new sequence 
        if (userParameters->getNoWeights() == false) {
            for (int j = 0; j < numSeqs; j++){
                if (aligned[j + 1] != 0){
                    (*seqWeight)[j] = static_cast<int>(treeWeight[j] * (*distMat)(j + 1, iseq));
                }
            }

            // Normalise the weights, such that the sum of the weights = INT_SCALE_FACTOR

            int sum = 0;
            for (int j = 0; j < numSeqs; j++) {
                if (aligned[j + 1] != 0){
                    sum += (*seqWeight)[j];
                }
            }

            if (sum == 0) {
                int j=0;
                for (j = 0; j < numSeqs; j++) {
                    (*seqWeight)[j] = 1;
                }

                sum = j;
            }
            for (int j = 0; j < numSeqs; j++) {
                if (aligned[j + 1] != 0) {
                    (*seqWeight)[j] = ((*seqWeight)[j] * INT_SCALE_FACTOR) / sum;
                    if ((*seqWeight)[j] < 1) {
                        (*seqWeight)[j] = 1;
                    }
                }
            }
        }

        int entries = 0;
        for (int j = 1; j <= numSeqs; j++){
            if (aligned[j] != 0) {
                group[j] = 1;
                entries++;
            }
            else if (iseq == j){
                group[j] = 2;
                entries++;
            }
        }
        
        alnPtr->addSeqWeight(seqWeight);
        aligned[iseq] = 1;
        
        ProfileAlignAlgorithm* alignAlgorithm = new MyersMillerProfileAlign;
        int score = alignAlgorithm->profileAlign(alnPtr, distMat, &group, aligned);
        delete alignAlgorithm;

        if (userParameters->getOutputOrder() == INPUT) {
            ++ix;
            newOutputIndex[iseq - 1] = iseq;
        }
        else {
            newOutputIndex[ix] = iseq;
            ++ix;
        }
    }
        
    alnPtr->addOutputIndex(&newOutputIndex);
    
    if(userParameters->getDisplayInfo()) {
        int alignmentScore = alnPtr->alignScore();
    }
    
    alnPtr->calculateMaxLengths(); // Mark change 20-6-07
    delete [] aligned;
    delete [] maxid;
    return (numSeqs);
}


