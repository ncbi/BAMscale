/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BAMstructs.h
 * Author: pongorls
 *
 * Created on November 28, 2018, 11:58 AM
 */

#ifndef BAMSTRUCTS_H
#define BAMSTRUCTS_H
#include "Definitions.h"

#ifdef __cplusplus
extern "C" {
#endif
    int CheckIndexShortFile(char *fname);
    int CheckIndexFile(char *fname);
    void DestroyBAMstruct(BAMFILES *head);
    BAMFILES *AddBAMstruct(char *BAMname, BAMFILES *head);
    void PrintBAMstructs(BAMFILES *head);
#ifdef __cplusplus
}
#endif
#endif /* BAMSTRUCTS_H */
