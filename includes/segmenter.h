/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   segmenter.h
 * Author: pongorls
 *
 * Created on November 29, 2018, 1:27 PM
 */

#ifndef SEGMENTER_H
#define SEGMENTER_H
#include "Definitions.h"
#include "main.h"

#ifdef __cplusplus
extern "C" {
#endif
    void DestroySegments(SEGMENTS *head);
    SEGMENTS *createSegment(void);
    int compare_float (const void * a, const void * b);
    int64_t CalculateGenSize(CHROMOSOMES *head);
    int64_t CalculateNonZeroBins(CHROMOSOMES *head, int sampleid);
    float *ConcatenateGenome(CHROMOSOMES *head, int64_t gensize, int sampleid);
    void Segmenting(CHROMOSOMES *head, CMDINPUT *cmd, int sampleid, float upper, float median, float lower);
    void Quantiles(CHROMOSOMES *head, int sampleid, CMDINPUT* cmd);
#ifdef __cplusplus
}
#endif
#endif /* SEGMENTER_H */
