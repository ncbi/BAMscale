/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CHROMstruct.h
 * Author: pongorls
 *
 * Created on November 28, 2018, 12:34 PM
 */

#ifndef CHROMSTRUCT_H
#define CHROMSTRUCT_H

#include "Definitions.h"
#include "main.h"

#ifdef __cplusplus
extern "C" {
#endif
    uint32_t *GetChrLens(CHROMOSOMES *head, int no_of_chrs);
    int CountNumberOfChromosomes(CHROMOSOMES *head);
    char **GetChromosomeNames(CHROMOSOMES *head, int no_of_chrs);
    float CalculateGenomeSize(CHROMOSOMES *head);
    CHROMOSOMES *AddCHROMstruct(CHROMOSOMES *head, char *name, int length, int no_of_samples, int threadID);
    void DestroyCHROMstruct(CHROMOSOMES *head, int no_of_samples);
    CHROMOSOMES *ImportChromosomeDataFromBAM(char *bamfile, int no_of_samples, int threads);
    void PrintChromosomes(CHROMOSOMES *head, int no_of_samples);
    void PrintBlacklistedChromosomes(CHROMOSOMES *head, int no_of_samples);
    CHROMOSOMES *ComputeBins(CHROMOSOMES *head, int binSize);
    CHROMOSOMES *AllocateBins(CHROMOSOMES *head, int no_of_samples);
    CHROMOSOMES *BlacklistChromosome(CHROMOSOMES *head, char *name);
    void BlacklistChromosomeFiles(CHROMOSOMES *head, char *filename);
    void DestroyChromCovStruct(CHRCOV *head);
    CHRCOV *CreateChromCovStruct(char *name, int id, int nbins);
    void DestroyRatioStruct(RATIOS *ptr);
    RATIOS *CreateRatioStruct(char *s1, char *s2, int id1, int id2);
#ifdef __cplusplus
}
#endif
#endif /* CHROMSTRUCT_H */
