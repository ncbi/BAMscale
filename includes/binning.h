/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   binning.h
 * Author: pongorls
 *
 * Created on November 28, 2018, 5:09 PM
 */

#ifndef BINNING_H
#define BINNING_H
#include "Definitions.h"
#include "main.h"

#ifdef __cplusplus
extern "C" {
#endif
    float *QuicksmoothenBins(float *carray, int smoothBins, int numBins);
    float *smoothenBins(float **carray, int smoothBins, int numBins);
    float *AddPseudoToZeroCov(float *coverage, int len);
    float *BinCoverage(int *coverage, int chr_len, int binSize, int nbins);
#ifdef __cplusplus
}
#endif
#endif /* BINNING_H */
