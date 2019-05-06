/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Inputs.h
 * Author: pongorls
 *
 * Created on December 11, 2018, 8:00 AM
 */

#ifndef INPUTS_H
#define INPUTS_H

#include "Definitions.h"
#include "main.h"

#define INPUTS_BASE "base"
#define INPUTS_GENOME "genome"
#define INPUTS_SCALED "scaled"
#define INPUTS_SCALE "scale"
#define INPUTS_UNSCALED "unscaled"
#define INPUTS_LOG2 "log2"
#define INPUTS_RFD "rfd"
#define INPUTS_END "endseq"
#define INPUTS_ENDR "endseqr"
#define INPUTS_REP "reptime"
#define INPUTS_RATIO "ratio"
#define INPUTS_SUBSTRACT "subtract"
#define INPUTS_READS "reads"
#define INPUTS_SINGLE "single"
#define INPUTS_PAIRED "paired"
#define INPUTS_AUTO "auto"
#define INPUTS_NO "no"
#define INPUTS_SMALLEST "smallest"
#define INPUTS_COV "cov"

#ifdef __cplusplus
extern "C" {
#endif
    CMDINPUT *CreateCMDinput(void);
    void PrintScaleMessage(char *pname);
    CMDINPUT *ScaleParser(int argc, char **argv);
    void PrintMultiCovMessage(char *pname);
    CMDINPUT *MultiCovParser(int argc, char **argv);
    void DestroyCMDinput(CMDINPUT *ptr);
#ifdef __cplusplus
}
#endif

#endif /* INPUTS_H */
