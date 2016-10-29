/**
 * Author: Andreas Wilm
 *
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.
 */
#ifndef STATS_H
#define STATS_H


#include <string>
#include <stdio.h>
#include "../alignment/Alignment.h"


namespace clustalw
{
    
using namespace std;


class Stats {
    
public:
    /* Functions */
    Stats();
    ~Stats();

    void setStatsFile(string f) {logfilename=f;};
    string getStatsFile() {return logfilename;};

    void setEnabled(bool b) {enabled=b;};
    bool isEnabled() {return enabled;};

    void logCmdLine(int argc, char **argv);

    void logInputSeqStats(Alignment *alnObj);// a bit confusing that this is
                                             // also an alignment...
    void logAlignedSeqStats(Alignment *alnObj);

    /* Attributes */
    
private:
    /* Functions */
#ifdef HAVE_MHASH_H
    char * Md5Hash(const char *thread);
    string Md5ForSeq(Alignment *alnObj, int s);
    string ConcatInputHash(Alignment *alnObj);
#endif
    /* adopted from Sean Eddy'ssquid:aligneval.c */
    float pairwiseIdentity(Alignment *alnObj, int s1, int s2);


    /* Attributes */
    string logfilename;
    bool enabled;
};
}
#endif
