/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.h
 * Author: pongorls
 *
 * Created on November 28, 2018, 1:21 PM
 */

#ifndef MAIN_H
#define MAIN_H
#include "Definitions.h"

#ifdef __cplusplus
extern "C" {
#endif
    void FreeAllocatedData(void);
    void ComputeCoverageGenome(char *infile, char *infile2, int binSize);
    void ComputeCoverageChIPpeak(CMDINPUT *cmd);
    void NormalizeBAMS(CMDINPUT *cmd);
    void PrintUsage(char *pname);
    int main();



#ifdef __cplusplus
}
#endif

#endif /* MAIN_H */

