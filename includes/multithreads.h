/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   multithreads.h
 * Author: pongorls
 *
 * Created on November 30, 2018, 12:04 PM
 */

#ifndef MULTITHREADS_H
#define MULTITHREADS_H
#include "Definitions.h"
#include "main.h"

#ifdef __cplusplus
extern "C" {
#endif
    void DestroyThreadStruct(THREADS **head, int no_of_threads);
    THREADS *CreateThreadStruct(char *chrname);
    THREADS *AddElement(THREADS *head, char *chrname);
    void * AnalyzeSample(void * voidA);
    void CalculateCoverageOfChromosomeBinsMultithreaded(CHROMOSOMES *head, BAMFILES *bhead, int paired_end, int bin_size, int pseudocount, int no_of_threads);
#ifdef __cplusplus
}
#endif
#endif /* MULTITHREADS_H */
