/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   binning.c
 * Author: pongorls
 * 
 * Created on November 28, 2018, 5:09 PM
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <htslib/sam.h>

#include "BAMcoverage.h"
#include "main.h"
#include "CHROMstruct.h"
#include "BAMstructs.h"
#include "binning.h"

float *QuicksmoothenBins(float *carray, int smoothBins, int numBins) {
    float *smoothed = (float *)calloc(numBins + 1, sizeof(float));
    int i = 0;
    int j = 0;
    int binstart = 0;
    int binend = 0;
    float binmean = 0;
    int runsum_state = 0;
    float runsum = -1;
    
    for(i = 0; i < numBins; i++) {
        binmean = 0;
        binstart = i - smoothBins;
        binend = i + smoothBins;
        runsum_state = 1;
        
        if(binstart < 1) {
            if(binstart < 0)
                binstart = 0;
            
            runsum_state = 0;
        }
        
        if(binend >= numBins) {
            binend = numBins - 1;
            runsum_state = 0;
        }
        
        if(binend - binstart > 0) {           
            if(runsum_state == 0) {
                for(j = binstart; j <= binend; j++) {
                    binmean += carray[j];
                }
            }
            
            else {
                binmean = runsum + carray[binend] - carray[binstart - 1];
            }
            
            runsum = binmean;
            
            if(binmean != 0) {
                binmean = binmean / (float)(binend - binstart);
            }
            
            smoothed[i] = binmean;
        }
    }
    
    if(carray)
        free(carray);
    
    return smoothed;
}


float *smoothenBins(float **carray, int smoothBins, int numBins) {
    float *smoothed = (float *) calloc(numBins, sizeof (float));
    int i = 0;
    int j = 0;
    int binstart = 0;
    int binend = 0;
    float binmean = 0;

    for (i = 0; i < numBins; i++) {
        binmean = 0;
        binstart = i - smoothBins;
        binend = i + smoothBins;

        if (binstart < 0)
            binstart = 0;

        if (binend >= numBins)
            binend = numBins - 1;

        if (binend - binstart > 0) {
            for (j = binstart; j <= binend; j++) {
                binmean += (*carray)[j];
            }

            if (binmean != 0) {
                binmean = binmean / (float) (binend - binstart);
            }

            smoothed[i] = binmean;
        } else {
            smoothed[i] = 0;
        }
    }

    if (*carray)
        free(*carray);

    return smoothed;
}

float *AddPseudoToZeroCov(float *coverage, int len) {
    int i = 0;

    for (i = 0; i < len; i++) {
        coverage[i] += 1;
    }

    return coverage;
}

float *BinCoverage(int *coverage, int chr_len, int binSize, int nbins) {
    float *bincov = (float *) calloc(nbins+1, sizeof (float));
    int i = 0;
    int j = 0;
    int binstart = 0;
    int binend = 0;
    
    //for (i = 0; i < nbins; i++) {
    //    bincov[i] = 0;
    //}

    for (i = 0; i < nbins; i++) {
        if (i == 0) {
            binstart = 0;
            binend = binSize;
        } else {
            binstart = i*binSize;
            binend = (i + 1) * binSize;
        }

        if (binstart >= chr_len) {
            binstart = chr_len-1;
        }

        if (binend >= chr_len) {
            binend = chr_len-1;
        }

        for (j = binstart; j < binend; j++) {
            bincov[i] += (float) coverage[j];
        }

        if(bincov[i] > 0)
            bincov[i] = bincov[i] / ((float) (binend - binstart));
    }

    //bincov = AddPseudoToZeroCov(bincov, nbins);

    return bincov;
}