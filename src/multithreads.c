/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   multithreads.c
 * Author: pongorls
 * 
 * Created on November 30, 2018, 12:04 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <pthread.h>
#include <htslib/sam.h>


#include "Definitions.h"
#include "main.h"
#include "scale.h"
#include "CHROMstruct.h"
#include "segmenter.h"
#include "BAMcoverage.h"
#include "multithreads.h"
#include "binning.h"

void DestroyThreadStruct(THREADS **head, int no_of_threads) {
    int i = 0;

    for (i = 0; i < no_of_threads; i++) {
        if ((*head)[i].chrname)
            free((*head)[i].chrname);
    }

    if (*head)
        free(*head);
}

THREADS *CreateThreadStruct(char *chrname) {
    THREADS *ptr = (THREADS *) malloc(sizeof (THREADS));

    ptr->chrname = strdup(chrname);

    ptr->sample_id = -1;
    ptr->paired_end = 0;
    ptr->scale = 1.00;
    ptr->binSize = 0;
    ptr->pseudocount = 1;
    ptr->sample = NULL;
    ptr->next = NULL;

    return ptr;
}

THREADS *AddElement(THREADS *head, char *chrname) {
    THREADS *curr = head;

    if (head == NULL) {
        head = CreateThreadStruct(chrname);
        curr = head;
    } else {
        while (curr->next != NULL) {
            curr = curr->next;
        }

        curr->next = CreateThreadStruct(chrname);
    }

    return head;
}

THREADS **AssignChrToThreads(CHROMOSOMES *head, int no_of_threads) {
    THREADS **ptr = (THREADS **) malloc(no_of_threads * sizeof (THREADS *));
    CHROMOSOMES *curr = head;
    int i = 0;

    while (curr != NULL) {
        if (curr->blacklist == 0) {
            ptr[i] = AddElement(ptr[i], curr->name);
            i++;

            if (i >= no_of_threads) {
                i = 0;
            }
        }

        curr = curr->next;
    }

    /*for(i = 0; i < no_of_threads; i++) {
        tmp = ptr[i];
        
        while(tmp != NULL) {
            //printf("[ %d ]\t%s\n", i, tmp->chrname);
            tmp = tmp->next;
        }   
    }*/

    return ptr;
}

CHROMOSOMES *FindChrStruct(CHROMOSOMES * head, char *chrname) {
    while (head != NULL) {
        if (strcmp(chrname, head->name) == 0) {
            return head;
        }

        head = head->next;
    }

    return NULL;
}

void * AnalyzeSample(void * voidA) {
    /*THREADS *ptr = (THREADS *) voidA;
    THREADS *ptr_head = ptr;
    CHROMOSOMES *tmp = ptr_head->head;
    char *filename = ptr_head->sample;
    
    int *cov = NULL;
    
    samFile *fp_in = hts_open(filename,"r"); //open bam file
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
    hts_idx_t *idx = sam_index_load(fp_in,  filename);
    hts_itr_t *iter  = NULL;
    
    while(ptr != NULL) {
        tmp = FindChrStruct(ptr_head->head, ptr->chrname);
        if(tmp != NULL) {
            iter = bam_itr_querys(idx, bamHdr, tmp->name);
            printf("%s\n", tmp->name);
            cov = CalculateCoverage(fp_in, iter, aln, tmp->length, tmp->name, ptr->cmd);
            tmp->coverages[ptr_head->sample_id] = BinCoverage(cov, tmp->length, ptr_head->binSize, tmp->numberOfBins);
            tmp->coverages[ptr_head->sample_id] = scaleBins(tmp->coverages[ptr_head->sample_id], ptr_head->scale, tmp->numberOfBins, (float)ptr_head->pseudocount);
                
            if(cov)
                free(cov);
                
            cov = NULL;
        }
        
        ptr = ptr->next;
    }
    
    bam_destroy1(aln);
    sam_close(fp_in);
     */
    return NULL;
}

void CalculateCoverageOfChromosomeBinsMultithreaded(CHROMOSOMES *head, BAMFILES *bhead, int paired_end, int bin_size, int pseudocount, int no_of_threads) {
    THREADS **threadStruct = AssignChrToThreads(head, no_of_threads);
    BAMFILES *bcurr = bhead;
    pthread_t thread_id[no_of_threads];
    int i = 0;

    while (bcurr != NULL) {
        for (i = 0; i < no_of_threads; i++) {
            if (threadStruct[i]) {
                threadStruct[i]->sample = bcurr->name;
                threadStruct[i]->sample_id = bcurr->id;
                threadStruct[i]->chr = head;
                threadStruct[i]->paired_end = paired_end;
                threadStruct[i]->scale = bcurr->scale;
                threadStruct[i]->binSize = bin_size;
                threadStruct[i]->pseudocount = pseudocount;
            }
        }

        for (i = 0; i < no_of_threads; i++) {
            if (threadStruct[i])
                pthread_create(&thread_id[i], NULL, &AnalyzeSample, threadStruct[i]);
        }

        for (i = 0; i < no_of_threads; i++) {
            if (threadStruct[i])
                pthread_join(thread_id[i], NULL);
        }

        bcurr = bcurr->next;
    }


    //for(i = 0; i < no_of_threads; i++) {
    //    DestroyThreadStruct(threadStruct[i]);
    //}

    if (threadStruct)
        DestroyThreadStruct(threadStruct, no_of_threads);
}

RATIOS *CalculateRatiosAllMultithreaded(RATIOS *head, CHROMOSOMES *chead, BAMFILES *bhead, int no_of_samples, int min_per_bin_cov, int smoothbin, int binSize, char *chromsizes) {
    RATIOS *curr = head;
    CHROMOSOMES *ccurr = chead;
    BAMFILES *bcurr = bhead;
    CHRCOV *ptr;

    if (no_of_samples < 2)
        return NULL;

    bcurr = bhead->next;

    while (bcurr != NULL) {
        printf("\nComparing samples:\n\t%s\n\t%s\n\n", bhead->name, bcurr->name);
        ccurr = chead;

        if (curr == NULL)
            curr = CreateRatioStruct(bhead->name, bcurr->name, bhead->id, bcurr->id);

        else {
            curr->next = CreateRatioStruct(bhead->name, bcurr->name, bhead->id, bcurr->id);
            curr = curr->next;
        }

        while (ccurr != NULL) {
            if (ccurr->blacklist == 0) {
                curr->chrcovs = CalculateChromosomeRatio(ccurr, curr->chrcovs, bhead->id, bcurr->id, 1, min_per_bin_cov);
            }

            ccurr = ccurr->next;
        }

        if (curr->chrcovs != NULL) {
            ptr = curr->chrcovs;

            while (ptr != NULL) {
                ptr->ratio = smoothenBins(&ptr->ratio, smoothbin, ptr->nbins);
                ptr = ptr->next;
            }
        }

        if (chromsizes)
            PrintBigWigOrdered(curr, binSize, chromsizes);
            //PrintBedgraphOrdered(curr, binSize, chromsizes);

        else
            PrintBedgraph(curr, binSize);
        //Quantiles(curr->chrcovs);

        bcurr = bcurr->next;
    }

    return head;
}