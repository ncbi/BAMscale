/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CHROMstruct.c
 * Author: pongorls
 * 
 * Created on November 28, 2018, 12:34 PM
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>

#include "Definitions.h"
#include "CHROMstruct.h"
#include "main.h"

uint32_t *GetChrLens(CHROMOSOMES *head, int no_of_chrs) {
    uint32_t *chrlens = (uint32_t *) malloc(sizeof (uint32_t) * no_of_chrs);
    CHROMOSOMES *curr = head;
    int i = 0;

    while (curr != NULL) {
        chrlens[i] = (uint32_t) curr->length;
        i++;
        curr = curr->next;
    }

    return chrlens;
}

int CountNumberOfChromosomes(CHROMOSOMES *head) {
    CHROMOSOMES *curr = head;
    int i = 0;

    while (curr != NULL) {
        i++;
        curr = curr->next;
    }

    return i;
}

char **GetChromosomeNames(CHROMOSOMES *head, int no_of_chrs) {
    char **chrnames = NULL;
    int i = 0;
    CHROMOSOMES *curr = head;

    if (no_of_chrs < 1)
        return NULL;

    chrnames = (char **) malloc(sizeof (char **) * no_of_chrs);

    while (curr != NULL) {
        if (i < no_of_chrs) {
            chrnames[i] = strdup(curr->name);
            i++;
        } else {
            fprintf(stderr, "WARNINGS: disregarding %s, there are more chromosomes than specified?\n", curr->name);
        }

        curr = curr->next;
    }

    return chrnames;
}

float CalculateGenomeSize(CHROMOSOMES *head) {
    float genome_size = 0;

    while (head != NULL) {
        if (head->blacklist == 0 && head->length > 0)
            genome_size += (float) head->length;

        head = head->next;
    }

    return genome_size;
}

/**
 * Creates new CHROMOSOME structure at end of HEAD.
 * @param <b>head</b> is the CHROMOSOME linked list head, <b>name</b> is the chr name
 * @return <b>HEAD</b> pointer to the CHROMOSOME structure head
 */
CHROMOSOMES *AddCHROMstruct(CHROMOSOMES *head, char *name, int length, int no_of_samples, int threadID) {
    CHROMOSOMES *ptr = (CHROMOSOMES *) malloc(sizeof (CHROMOSOMES));
    CHROMOSOMES *curr = head;

    ptr->blacklist = 0;
    ptr->length = length;
    ptr->numberOfBins = -1;
    ptr->id = -1;
    ptr->allocated = 0;

    ptr->name = NULL;
    ptr->next = NULL;
    ptr->coverages = NULL;
    ptr->tid = threadID;
    ptr->name = strdup(name);
    ptr->idxreads = (int *) malloc(no_of_samples * sizeof (int));
    
    for(int i =0 ; i < no_of_samples; i++){
        ptr->idxreads[i] = 0;
    }

    if (head == NULL) {
        ptr->id = 0;
        head = ptr;
    } else {
        while (curr->next != NULL) {
            curr = curr->next;
        }

        curr->next = ptr;
        ptr->id = curr->id + 1;
    }

    return head;
}

/*
 * 
 */
void DestroyCHROMstruct(CHROMOSOMES *head, int no_of_samples) {
    CHROMOSOMES *curr = head;
    int i = 0;

    while (head != NULL) {
        curr = head;
        head = head->next;
        curr->next = NULL;

        if (curr->coverages) {
            for (i = 0; i < no_of_samples; i++) {
                if (curr->coverages[i]){
                    free(curr->coverages[i]);
                }
            }
            free(curr->coverages);
        }

        if (curr->idxreads)
            free(curr->idxreads);

        if (curr->name)
            free(curr->name);

        if (curr)
            free(curr);
    }
}

CHROMOSOMES *ImportChromosomeDataFromBAM(char *bamfile, int no_of_samples, int threads) {
    CHROMOSOMES *head = NULL;
    samFile *fp_in = hts_open(bamfile, "r");
    bam_hdr_t *hdr = sam_hdr_read(fp_in);
    int i = 0;
    int j = 0;

    for (i = 0; i < hdr->n_targets; i++) {
        head = AddCHROMstruct(head, hdr->target_name[i], (int) hdr->target_len[i], no_of_samples, j);

        j++;

        if (j >= threads) {
            j = 0;
        }
    }
    
    bam_hdr_destroy(hdr);
    sam_close(fp_in);
    return head;
}

void PrintChromosomes(CHROMOSOMES *head, int no_of_samples) {
    int i;
    CHROMOSOMES *curr = head;

    while (curr != NULL) {
        printf("Name: %s\n", curr->name);
        printf("\tID: %d\n", curr->id);
        printf("\tLength: %d\n", curr->length);
        printf("\tBlacklisted: %d\n", curr->blacklist);
        printf("\tNumber of bins: %d\n", curr->numberOfBins);
        printf("\tAllocated sample bins: %d\n", curr->allocated);
        printf("\tThread ID: %d\n", curr->tid);

        if (curr->idxreads) {
            for (i = 0; i < no_of_samples; i++)
                printf("\t\t[ %d ] = %d\n", i, curr->idxreads[i]);
        }
        curr = curr->next;
    }
}

void PrintBlacklistedChromosomes(CHROMOSOMES *head, int no_of_samples) {
    int i;
    CHROMOSOMES *curr = head;

    while (curr != NULL) {
        if (curr->blacklist == 1) {
            printf("Name: %s\n", curr->name);
            printf("\tID: %d\n", curr->id);
            printf("\tLength: %d\n", curr->length);
            printf("\tBlacklisted: %d\n", curr->blacklist);
            printf("\tNumber of bins: %d\n", curr->numberOfBins);
            printf("\tAllocated sample bins: %d\n", curr->allocated);

            if (curr->idxreads) {
                for (i = 0; i < no_of_samples; i++)
                    printf("\t\t[ %d ] = %d\n", i, curr->idxreads[i]);
            }
        }
        curr = curr->next;
    }
}

CHROMOSOMES *ComputeBins(CHROMOSOMES *head, int binSize) {
    CHROMOSOMES *curr = head;

    while (curr != NULL) {
        curr->numberOfBins = curr->length / binSize;
        curr = curr->next;
    }

    return head;
}

CHROMOSOMES *AllocateBins(CHROMOSOMES *head, int no_of_samples) {
    CHROMOSOMES *curr = head;
    int i = 0;

    if (no_of_samples <= 0) {
        printf("ERROR: no samples were specified??");
        FreeAllocatedData();
        exit(0);
    }

    while (curr != NULL) {
        if (curr->numberOfBins > -1 && curr->blacklist == 0) {
            curr->coverages = (float **) calloc(no_of_samples+1, sizeof (float **));

            if (curr->coverages == NULL) {
                printf("ERROR: could not allocate memory for bins at chr: %s\n", curr->name);
                FreeAllocatedData();
                exit(0);
            }

            for (i = 0; i < no_of_samples; i++) {
                /*curr->coverages[i] = (float *) malloc(sizeof (float) * no_of_samples);

                if (curr->coverages[i] == NULL) {
                    printf("ERROR: could not allocate memory for bins at chr: %s for sample no. [ %d ]\n", curr->name, i);
                    FreeAllocatedData();
                    exit(0);
                }*/

                curr->allocated++;
            }
        } else {
            curr->coverages = NULL;
        }

        curr = curr->next;
    }

    return head;
}

CHROMOSOMES *BlacklistChromosome(CHROMOSOMES *head, char *name) {
    CHROMOSOMES *curr = head;
    int found = 0;

    while (curr != NULL) {
        if (strcmp(name, curr->name) == 0) {
            curr->blacklist = 1;
            curr->tid = -1;
            found++;
        }

        curr = curr->next;
    }

    if (found == 0) {
        printf("WARNING: \"%s\" chromosome not found and could not be blacklisted\n", name);
    }

    return head;
}

void BlacklistChromosomeFiles(CHROMOSOMES *head, char *filename) {
    FILE *handler = fopen(filename, "r");
    char line[BUFSIZ];
    char *pos;

    while (fgets(line, sizeof (line), handler)) {
        if ((pos = strchr(line, '\n')) != NULL)
            *pos = '\0';

        head = BlacklistChromosome(head, line);
    }

    fclose(handler);
    //return head;
}

void DestroyChromCovStruct(CHRCOV *head) {
    CHRCOV *curr = head;

    while (head != NULL) {
        curr = head;
        head = head->next;

        if (curr->name)
            free(curr->name);

        if (curr->ratio)
            free(curr->ratio);

        if (curr)
            free(curr);
    }
}

CHRCOV *CreateChromCovStruct(char *name, int id, int nbins) {
    CHRCOV *ptr = (CHRCOV *) malloc(sizeof (CHRCOV));

    ptr->name = strdup(name);

    ptr->id = id;
    ptr->nbins = nbins;

    ptr->next = NULL;
    ptr->ratio = NULL;

    /*if(ptr->ratio == NULL) {
        printf("ERROR: could not allocate memory for bins at chr: %s for sample no. [ %d ]\n", ptr->name, id);
        FreeAllocatedData();
        exit(0); 
    }*/

    return ptr;
}

void DestroyRatioStruct(RATIOS *ptr) {
    if (ptr->sample1)
        free(ptr->sample1);

    if (ptr->sample2)
        free(ptr->sample2);

    if (ptr)
        free(ptr);

    if (ptr->chrcovs)
        DestroyChromCovStruct(ptr->chrcovs);

    ptr = NULL;
}

RATIOS *CreateRatioStruct(char *s1, char *s2, int id1, int id2) {
    RATIOS *ptr = (RATIOS *) malloc(sizeof (RATIOS));
    ptr->chrcovs = NULL;
    ptr->s1 = id1;
    ptr->s2 = id2;

    ptr->sample1 = strdup(s1);

    ptr->sample2 = strdup(s2);

    return ptr;
}