/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   scale.h
 * Author: pongorls
 *
 * Created on November 28, 2018, 3:54 PM
 */

#ifndef SCALE_H
#define SCALE_H
#include "Definitions.h"
#include "main.h"

#ifdef __cplusplus
extern "C" {
#endif
    float *scaleBins(float *carray, float scale, int nbins, float pseudocount);
    void ScaleToSmallest(BAMFILES *head);
    void NoScale(BAMFILES *head);
    void ScaleToGenomeSize(BAMFILES *head, CHROMOSOMES *chead);
    void ScaleGenomeCoverage(BAMFILES *head, CHROMOSOMES *chead);
    BAMFILES *ComputeSamplescales(BAMFILES *head, CHROMOSOMES *chead, int scale);
    float *logTwoCoverageRatio(float *cov1, float *cov2, int nbins, float min_per_bin_cov);
    float *OKseqRFD(float *cov1, float *cov2, int nbins, float min_per_bin_cov);
    float *SubtractCoverage(float *cov1, float *cov2, int nbins, float min_per_bin_cov);
    float *CoverageRatio(float *cov1, float *cov2, int nbins, float min_per_bin_cov);
    float *SignedCoverageRatio(float *cov1, float *cov2, int nbins, int min_per_bin_cov);
    CHRCOV *CalculateChromosomeRatio(CHROMOSOMES *curr, CHRCOV *chead, int s1, int s2, int ratioType, int min_per_bin_cov);
    void PrintBedgraph(RATIOS *ptr, int binSize);
    char *returnChrName(char *input);
    void PrintBedgraphOrdered(RATIOS *ptr, int binSize, char *chromfile);
    void PrintBigWigOrdered(RATIOS *ptr, int binSize, char *chromfile);
    RATIOS *CalculateRatiosAll(RATIOS *head, CHROMOSOMES *chead, BAMFILES *bhead, int no_of_samples, int min_per_bin_cov, int smoothbin, int binSize, char *chromsizes);
#ifdef __cplusplus
}
#endif
#endif /* SCALE_H */
