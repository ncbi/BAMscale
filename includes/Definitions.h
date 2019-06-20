/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Definitions.h
 * Author: pongorls
 *
 * Created on November 28, 2018, 11:56 AM
 */

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
    typedef struct segment {
        int start;
        int end;
        int value;
        
        struct segment *next;
        struct segment *prev;
    } SEGMENTS;
    
    
    typedef struct peak {
        char *coord;
        int id;
        char *chr;
        char *start_str;
        char *end_str;
        int strand;
        
        int start;
        int end;
        int length;
        
        int nbins;
        int binSize;
        int noSamples;
        
        int tid;
        
        int *read_cov;
        float *normalized;
        float **cov;
        struct peak *next;
    } PEAK;
    
    typedef struct chrcov {
        char *name;
        int id;
        int nbins;
        float *ratio;
        
        struct chrcov *next;
    } CHRCOV;
    
    typedef struct ratios {
        char *sample1;
        char *sample2;
        int s1;
        int s2;
        
        CHRCOV *chrcovs;
        struct ratios *next;
    } RATIOS;
    
    typedef struct BEDcoords {
        char *coord;
        char *chr;
        char *start_str;
        char *end_str;
        float *normcov;
        
        int start;
        int end;
        int summit;
        int strand;
        struct BEDcoords *next;
    } BEDCOORDS;
    
    typedef struct BAMfiles{
        char *name;
        char *shortname;
        int id;
        int read_coverage;
        int filtered_reads;
        double base_coverage;
        float scale;
        float genome_scale;
        struct BAMfiles *next;
    } BAMFILES;
    
    typedef struct chromosomes {
        char *name;
        int id;
        int length;
        int accept;
        int blacklist;
        float  **coverages;
        int *idxreads;
        int numberOfBins;
        int allocated;
        int nonzerobins;
        int tid;
        struct chromosomes *next;
    } CHROMOSOMES;

    typedef struct cmdinput {
        char *bedfile;
        int no_of_samples;
        BAMFILES *bamfiles;
        CHROMOSOMES *chr;
        int mapq;
        int removeduplicates;
        int nounproper;
        int remove_unmapped_pair;
        int fragment_count_mode;
        int fragment_size_filter;
        int fragment_size;
        int strand;
        int libtype;
        int threads;
        int min_insert_size;
        int max_insert_size;
        char *outdir;
        char *outprefix;
        char *blacklist_file;
        char *blacklist_bed;
        char **argv;
        int argc;
        int genome_coverage;
        int strandsplit;
        char *custom_scale;
        
        char *normtype;
        char *scale;
        char *operation;
        
        int filtDiffChr;
        int filtInsSize;
        
        int binSize;
        int binSizeChange;
        int smoothBin;
        int smoothBinChange;
        int tracksmooth;
    } CMDINPUT;
    
    typedef struct threads {
        int pid;
        char *chrname;
        char *sample;
        int sample_id;
        int paired_end;
        float scale;
        int binSize;
        int pseudocount;
        int strand;
        CHROMOSOMES *chr;
        PEAK *phead;
        CMDINPUT *cmd;
        BAMFILES *bamfile;
        struct threads *next;
    } THREADS;
#ifdef __cplusplus
}
#endif

#endif /* DEFINITIONS_H */

