/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "RootedTreeOutput.h"

namespace clustalw
{

RootedTreeOutput::RootedTreeOutput(SeqInfo* seqInfo)
{
    firstSeq = seqInfo->firstSeq;
    lastSeq = seqInfo->lastSeq;
    numSeqs = seqInfo->numSeqs;  
}

void RootedTreeOutput::printPhylipTree(RootedGuideTree* tree, ofstream* ptrToFile, Alignment *alignPtr,
                                    DistMatrix* distMat)
{
    if(!ptrToFile || !ptrToFile->is_open())
    {
        return;
    }
    
    // If we have only 2 sequences, use the distances in the distMat
    if (lastSeq - firstSeq + 1 == 2)
    {
        (*ptrToFile) << "(" << alignPtr->getName(firstSeq) << ":" << fixed << setprecision(5) 
        << (*distMat)(firstSeq, firstSeq + 1) << "," << alignPtr->getName(firstSeq + 1) 
        << ":" << fixed << setprecision(5)  << (*distMat)(firstSeq, firstSeq + 1);
        return ;
    }
    
    (*ptrToFile) << "(\n";     
    phylipTraverse(ptrToFile, alignPtr, tree->getRoot());
    (*ptrToFile) << ");\n";
}

void RootedTreeOutput::printNexusTree(RootedGuideTree* tree, ofstream* ptrToFile, 
                                      Alignment *alignPtr, DistMatrix* distMat)
{
    if(!ptrToFile || !ptrToFile->is_open())
    {
        return;
    }
        
    (*ptrToFile) << "#NEXUS\n\n";

    (*ptrToFile) << "BEGIN TREES;\n\n";
    (*ptrToFile) << "\tTRANSLATE\n";
            
    for(int j = 1; j < numSeqs; j++) 
    {
        (*ptrToFile) << "\t\t" << j << "\t" << alignPtr->getName(j) <<",\n";
    }
    (*ptrToFile) << "\t\t" << numSeqs << "\t" << alignPtr->getName(numSeqs) << "\n";
    (*ptrToFile) << "\t\t;\n";

    (*ptrToFile) << "\tUTREE PAUP_1= ";
    
    // IF we have only 2 seqs
    if (lastSeq - firstSeq + 1 == 2)
    {
        (*ptrToFile) << "(" << alignPtr->getName(firstSeq) << ":" << fixed << setprecision(5) 
        << (*distMat)(firstSeq, firstSeq + 1) << "," << alignPtr->getName(firstSeq + 1) 
        << ":" << fixed << setprecision(5)  << (*distMat)(firstSeq, firstSeq + 1);
    }
    else
    {                        
        (*ptrToFile) << "(";
        nexusTraverse(ptrToFile, alignPtr, tree->getRoot());
    }
    (*ptrToFile) << ");\n";
    (*ptrToFile) << "\nENDBLOCK;\n";
}

/**
 * PRIVATE FUNCTIONS
 */ 
 
void RootedTreeOutput::phylipTraverse(ofstream* ptrToFile, Alignment *alignPtr, Node* t)
{
    if(!ptrToFile)
    {
        return;
    }    
    if(t != 0)
    {
        if(t->isLeafNode())
        {
            if(alignPtr)
            {
                (*ptrToFile) << alignPtr->getName(t->getSeqNum()) << ":" << t->getHeight();
            }
            else
            {
                (*ptrToFile) << t->getSeqNum() << ":" << t->getHeight();
            }
        }
        else // Internal node
        {
            (*ptrToFile) << "(\n";
            phylipTraverse(ptrToFile, alignPtr, t->getLeft());
            (*ptrToFile) << ",\n";
            phylipTraverse(ptrToFile, alignPtr, t->getRight());
            (*ptrToFile) << "):" << t->getHeight();
        }
    }
}

void RootedTreeOutput::nexusTraverse(ofstream* ptrToFile, Alignment *alignPtr, Node* t)
{
    if(t != 0)
    {
        if(!t->isLeafNode()) // Internal node
        {
            (*ptrToFile) << "(";
            nexusTraverse(ptrToFile, alignPtr, t->getLeft());
            (*ptrToFile) << ",";
            nexusTraverse(ptrToFile, alignPtr, t->getRight());
            (*ptrToFile) << "):" << t->getHeight();                    
        }
        else // Leaf node
        {
            if(alignPtr)
            {
                (*ptrToFile) << alignPtr->getName(t->getSeqNum()) << ":" << t->getHeight();
            }
            else
            {
                (*ptrToFile) << t->getSeqNum() << ":" << t->getHeight();
            }
        }
    }

}

}
