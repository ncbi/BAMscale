/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Writer.h
 * Author: pongorls
 *
 * Created on December 19, 2018, 2:27 PM
 */

#ifndef WRITER_H
#define WRITER_H

#include "Definitions.h"
#include "main.h"

#ifdef __cplusplus
extern "C" {
#endif
    char *returnRNAfilename(CMDINPUT *cmd);
    void PrintScaledBigWig(CMDINPUT *cmd, BAMFILES *curr, char *sfile);
#ifdef __cplusplus
}
#endif

#endif /* WRITER_H */
