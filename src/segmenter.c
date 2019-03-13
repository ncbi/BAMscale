/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   segmenter.c
 * Author: pongorls
 * 
 * Created on November 29, 2018, 1:27 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "Definitions.h"
#include "main.h"
#include "scale.h"
#include "CHROMstruct.h"
#include "segmenter.h"

void DestroySegments(SEGMENTS *head) {
    SEGMENTS *curr = head;
    
    while(head != NULL) {
        curr = head;
        head = head->next;
        
        free(curr);
    }
}
SEGMENTS *createSegment(void) {
    SEGMENTS *ptr = (SEGMENTS *) malloc(sizeof(SEGMENTS));
    ptr->next = NULL;
    ptr->prev = NULL;
    ptr->start = -1;
    ptr->end = -1;
  
    return ptr;
}

int compare_float (const void * a, const void * b) {
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}


int64_t CalculateGenSize(CHROMOSOMES *head) {
    int64_t gensize = 0;
    CHROMOSOMES *curr = head;

    while(curr != NULL) {
        if(curr->blacklist == 0) {
            gensize += (int64_t)curr->numberOfBins;      
        }
        
        curr = curr->next;
    }
    return gensize;
}

int64_t CalculateNonZeroBins(CHROMOSOMES *head, int sampleid) {
    int64_t nzbins = 0;
    CHROMOSOMES *curr = head;
    int i;
    
    while(curr != NULL) {
        if(curr->blacklist == 0) {
            curr->nonzerobins = 0;
            
            for(i = 0; i < curr->numberOfBins; i++) {
                if(curr->coverages[sampleid][i] != 0) {
                    curr->nonzerobins++;
                }
            }     
        }
        
        curr = curr->next;
    }
    
    curr = head;
    
    while(curr != NULL) {
        if(curr->blacklist == 0) {
            nzbins += curr->nonzerobins;
        }
        
        curr = curr->next;
    }
    
    return nzbins;
}

float *ConcatenateGenome(CHROMOSOMES *head, int64_t gensize, int sampleid) {
    CHROMOSOMES *curr = head;
    int64_t i = 0;
    int64_t j = 0;
    float *genbins = NULL;
    
    
    if(gensize < 1)
        return NULL;
    
    genbins = (float *) malloc(gensize * sizeof(float));
    
    if(genbins == NULL) {
        printf("ERROR: could not allocate memory for genome (at quantile calculation)\n");
        FreeAllocatedData();
        exit(0); 
    }
        
    while(curr != NULL) {
        if(curr->blacklist == 0) {
            
            for(i = 0; i < curr->numberOfBins; i++) {
                if(curr->coverages[sampleid][i] != 0) {
                   genbins[j] = curr->coverages[sampleid][i];
                   j++; 
                }
            }
        }
        
        curr = curr->next;
    }
    
    
    
    return genbins;
}

void Segmenting(CHROMOSOMES *head, CMDINPUT *cmd, int sampleid, float upper, float median, float lower) {
    SEGMENTS *segment = NULL;
    SEGMENTS *segmenthead = NULL;
    SEGMENTS *tmp = NULL;
    SEGMENTS *rtmp = NULL;
    CHROMOSOMES *curr = head;
    int prevstate = -1;
    int currstate = -1;
    int start;
    int end;
    int i = 0;
    FILE *fp_s = fopen("strong.bed", "w+");
    FILE *fp_sm = fopen("med_strong.bed", "w+");
    FILE *fp_wm = fopen("med_weak.bed", "w+");
    FILE *fp_w = fopen("weak.bed", "w+");
    int minsize = 200;
    printf("Minsize: %d\n", minsize);
    
    while(curr != NULL) {
        if(curr->blacklist == 0) {
            currstate = -1;
            prevstate = -1;
            start = 0;
            end = 0;
            segment = NULL;

            for(i = 0; i < curr->numberOfBins; i++) {
                //printf("%d: %f\t%f\n", i, curr->coverages[sampleid][i], median);
                
                if(curr->coverages[sampleid][i] == 0) {
                    currstate = 0;
                }
                
                else if(curr->coverages[sampleid][i] > upper) {
                   currstate = 4; 
                }
                
                else if(curr->coverages[sampleid][i] > median) {
                    currstate = 3;
                }
                
                else if(curr->coverages[sampleid][i] > lower) {
                    currstate = 2;
                }
                
                else {
                   currstate = 1; 
                }
                
                if(prevstate == -1) {
                    start = 0;
                    prevstate = currstate;
                }
                
                if(currstate == prevstate) {
                    end = i;
                }

                else {
                    if(segmenthead == NULL) {
                        segmenthead = createSegment();
                        segment = segmenthead;
                    }
                    
                    else {
                        segment->next = createSegment();
                        segment->next->prev = segment;
                        
                        segment = segment->next;
                    }
                    
                    segment->value = prevstate;
                    segment->start = start;
                    segment->end = end;
                    
                    prevstate = currstate;
                    start = i-1;
                    end = i;
                }
            }
            
            segment = segmenthead;
            
            while(segment != NULL) {
                if(segment->end - segment->start >= minsize) {
                    tmp = segment->next;
                    
                    while(tmp != NULL && tmp->end - tmp->start < minsize) {
                        rtmp = tmp;
                        tmp = tmp->next;
                        free(rtmp);
                    }
                    
                    if(tmp == NULL) {
                        segment->end = curr->numberOfBins;
                        segment->next = NULL;
                    }
                    
                    else {
                        if(tmp != segment->next) {
                            /*if(tmp->start - segment->end < minsize) {
                                if(tmp->value == segment->value) {
                                    segment->end = tmp->start;
                                }
                            }*/
                            
                            
                            segment->next = tmp;
                            tmp->prev = segment;
                        }
                    }
                }
                
                segment = segment->next;
                
            }
            
            segment = segmenthead;
            
            while(segment != NULL) {
                if(segment->next) {
                    if(segment->end != segment->next->start) {
                        if(segment->end - segment->start > segment->next->end - segment->next->start)
                            segment->end = segment->next->start;
                        
                        else
                            segment->next->start = segment->end;
                    }
                }
                
                segment = segment->next;
            }
            
            segment = segmenthead;
            
            while(segment != NULL) {
                //if(segment->end - segment->start > minsize) {
                    if(segment->value == 4)
                        fprintf(fp_s, "%s\t%d\t%d\n", curr->name, segment->start*cmd->binSize, segment->end*cmd->binSize);

                    if(segment->value == 3)
                        fprintf(fp_sm, "%s\t%d\t%d\n", curr->name, segment->start*cmd->binSize, segment->end*cmd->binSize);

                    if(segment->value == 2)
                        fprintf(fp_wm, "%s\t%d\t%d\n", curr->name, segment->start*cmd->binSize, segment->end*cmd->binSize);

                    if(segment->value == 1)
                        fprintf(fp_w, "%s\t%d\t%d\n", curr->name, segment->start*cmd->binSize, segment->end*cmd->binSize);
                //}
                    //if(segment->end - segment->start >= 50 && segment->value > 2)
                    //printf("%d\t%d\t%d\n", segment->start*5, segment->end*5, segment->value);
                    segment = segment->next;
            }
            
            printf("Done\n");
            
            DestroySegments(segmenthead);
            segment = NULL;
            segmenthead = NULL;
        }
        
        curr = curr->next;
    }
    
    //DestroySegments(segment);
    
    fclose(fp_s);
    fclose(fp_sm);
    fclose(fp_wm);
    fclose(fp_w);
}

void Quantiles(CHROMOSOMES *head, int sampleid, CMDINPUT *cmd) {
    //int64_t gensize = CalculateGenSize(head);
    int64_t gensize = CalculateNonZeroBins(head, sampleid);
    fprintf(stderr, "Genome size: %" PRIu64 "\n", gensize);
 
    float *genbin = ConcatenateGenome(head, gensize, sampleid);
    
    float upper = -1;
    float median = -1;
    float lower = -1;
    
    printf("Sorting genome ( %d ) \n", sampleid);
    
    qsort(genbin, gensize, sizeof(float), compare_float);
    
    upper = genbin[(int)round(gensize*0.75)];
    median = genbin[(int)round(gensize*0.5)];
    lower = genbin[(int)round(gensize*0.25)];
            
    printf("Quantiles:\n");
    printf("\tUpper: %f\n", upper);
    printf("\tMedian: %f\n", median);
    printf("\tLower: %f\n", lower);

    Segmenting(head, cmd, sampleid, upper, median, lower);
    
    if(genbin)
        free(genbin);
}