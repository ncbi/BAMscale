/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BEDstruct.h
 * Author: pongorls
 *
 * Created on December 10, 2018, 8:02 AM
 */

#ifndef BEDSTRUCT_H
#define BEDSTRUCT_H

#include "Definitions.h"
#include "main.h"

#ifdef __cplusplus
extern "C" {
#endif
    char *BEDtoString(char *chr, int start, int end);
    int Read_filter_MultiCov(bam1_t *read, int paired_end);
    void DeleteBEDs(PEAK *head);
    PEAK *CreateBEDentry(void);
    void AllocateReadCovs(PEAK *head, int no_of_samples);
    void AllocateCovs(PEAK *head);
    PEAK *AddBEDentry(PEAK *curr, char *BEDentry, int tid);
    PEAK *ReadBED(char *BEDfilename, int nthreads);
    void GetBEDCoveragesBAM(BAMFILES *bhead, PEAK *beds, int paired_end);
    void *GetBEDFragmentCoveragesBAMmultithread(void * voidA);
    void *GetBEDCoveragesBAMmultithread(void * voidA);
    void MultiCoverage(BAMFILES *bhead, PEAK *head, CMDINPUT *cmd);
    void CalculateFPKM(BAMFILES *bhead, PEAK *head);
    void CalculateLibScaled(BAMFILES *bhead, PEAK *head);
    void CalculateTPM(BAMFILES *bhead, PEAK *head);
    void WriteMultiCovsRaw(BAMFILES *bhead, PEAK *head, int no_of_samples, char *outfile);
    void WriteMultiCovsNormalized(BAMFILES *bhead, PEAK *head, int no_of_samples, char *outfile);
#ifdef __cplusplus
}
#endif

#endif /* BEDSTRUCT_H */
