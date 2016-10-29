/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef RSFFILEPARSER_H
#define RSFFILEPARSER_H

#include <string>
#include "FileParser.h"

namespace clustalw
{

class RSFFileParser : public FileParser
{
    public:
        /* Functions */
        RSFFileParser(string filePath);
        virtual vector<Sequence> getSeqRange(int firstSeq, int num);
        virtual Sequence getSeq(int seqNum);
        virtual int countSeqs();
        virtual void getSecStructure(vector<char>& gapPenaltyMask, 
                                     vector<char>& secStructMask, string& secStructName, 
                                     int &structPenalties, int length); 

        /* Attributes */

    private:
        /* Functions */
        void getRSFFeature(char* line, vector<char>& secStructMask, int length);
        bool keyword(char *line,char *code);
        /* Attributes */
        string fileName;
};

}
#endif


