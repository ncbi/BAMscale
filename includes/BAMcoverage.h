/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BAMcoverage.h
 * Author: pongorls
 *
 * Created on November 28, 2018, 12:27 PM
 */

#ifndef BAMCOVERAGE_H
#define BAMCOVERAGE_H

#include "Definitions.h"
#include "main.h"
#include <htslib/sam.h>

#ifdef __cplusplus
extern "C" {
#endif
    int ReadStrand(bam1_t *read, int paired_end);
    int DetectLibraryType(BAMFILES *bhead);
    int Read_filter(bam1_t *read, CMDINPUT *cmd);
    CHROMOSOMES *AddIDXcoverage(char *name, int coverage, int id, CHROMOSOMES *head);
    void GetChromosomeCoveragesIDX(CHROMOSOMES *head, BAMFILES *bhead);
    void GetGenomeCoveragesIDX(CHROMOSOMES *head, BAMFILES *bhead);
    void CalculateCoverageOfReads(samFile *fp_in, hts_itr_t *iter, bam1_t *aln, int chrsize, char *chrname, CMDINPUT *cmd, BAMFILES *bamcurr);
    void *GetGenomeReadCoveragemultithread(void * voidA);
    void MultiGenomeReadCoverage(CMDINPUT *cmd, CHROMOSOMES *chr);
    void GetChromosomeCoveragesBAM(CHROMOSOMES *head, BAMFILES *bhead, CMDINPUT *cmd);
    char *BEDentryToCoord(char *input);
    char *BEDentryChr(char *input);
    void SubtractBlacklistedBEDS(char *filename, CHROMOSOMES *head, BAMFILES *bhead, int paired_end);
    int *CalculateCoverage(samFile *fp_in, hts_itr_t *iter, bam1_t *aln, int chrsize, char *chrname, CMDINPUT *cmd, BAMFILES *bamcurr);
    void GetGenomeCoverageRNA(CMDINPUT *cmd, CHROMOSOMES *head, char *outfile);
    void *GetGenomeCoveragemultithread(void * voidA);
    void MultiGenomeCoverage(CMDINPUT *cmd, CHROMOSOMES *chr);
    void *GetGenomeBaseCoveragemultithread(void * voidA);
    void MultiGenomeBaseCoverage(CMDINPUT *cmd, CHROMOSOMES *chr);
    void *ScaleBinsmultithread(void * voidA);
    void MultiGenomeScaler(CMDINPUT *cmd, CHROMOSOMES *chr);
    void *SmoothBinsmultithread(void * voidA);
    void MultiGenomeSmoother(CMDINPUT *cmd, CHROMOSOMES *chr);
    void *TransformBinsmultithread(void * voidA);
    void MultiGenomeTransform(CMDINPUT *cmd, CHROMOSOMES *chr);
    void CalculateCoverageOfChromosomeBins(CHROMOSOMES *head, BAMFILES *bhead, int paired_end, int bin_size, int pseudocount, CMDINPUT *cmd);
    void SmoothenAllChromosomeBins(CHROMOSOMES *head, BAMFILES *bhead, int smoothBinNum);
#ifdef __cplusplus
}
#endif
#endif /* BAMCOVERAGE_H */
