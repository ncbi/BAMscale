/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Writer.c
 * Author: pongorls
 * 
 * Created on December 19, 2018, 2:27 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libgen.h>

#include <bigWig.h>

#include "Writer.h"
#include "Definitions.h"
#include "main.h"
#include "scale.h"
#include "CHROMstruct.h"
#include "segmenter.h"
#include "binning.h"
#include "Inputs.h"

void PrintScaledBigWig(CMDINPUT *cmd, BAMFILES *curr, char *sfile) {
    char **chrnames = NULL;
    uint32_t *chrlens = NULL;
    uint32_t start = 0;
    int no_of_chrs = CountNumberOfChromosomes(cmd->chr);
    char *outfile = NULL;
    int fnamelen = 0;
    bigWigFile_t *fp = NULL;
    int i = 0;
    CHROMOSOMES *chr = cmd->chr;
    float *intervals = NULL;

    chrnames = GetChromosomeNames(cmd->chr, no_of_chrs);
    chrlens = GetChrLens(cmd->chr, no_of_chrs);

    if (cmd->outdir != NULL)
        fnamelen += strlen(cmd->outdir);

    fnamelen += strlen(curr->shortname);

    if(cmd->strandsplit == 1)
        fnamelen += strlen(".positive");
    
    if (sfile != NULL)
        fnamelen += strlen(sfile);

    fnamelen += 50;

    outfile = (char *) calloc((fnamelen*2 + 1), sizeof (char));

    if (cmd->outdir != NULL)
        strcpy(outfile, cmd->outdir);

    else
        strcpy(outfile, "./");

    strcat(outfile, "/");

    if (sfile == NULL) {
        strcat(outfile, curr->shortname);
        
        if(cmd->strandsplit == 1) {
            if(cmd->strand == 1)
                strcat(outfile, ".positive");
            
            if(cmd->strand == -1)
                strcat(outfile, ".negative");
        }
    }
    else {
        strcat(outfile, curr->shortname);
        
        if(strcmp(cmd->operation, INPUTS_END) != 0 && strcmp(cmd->operation, INPUTS_ENDR) != 0 && strcmp(cmd->operation, INPUTS_RFD) != 0) {
            strcat(outfile, "_vs_");
            strcat(outfile, sfile);
        }
        
        if(strcmp(cmd->operation, INPUTS_END) == 0 || strcmp(cmd->operation, INPUTS_ENDR) == 0) {
            strcat(outfile, ".log2");
        }
    }

    strcat(outfile, ".");
    strcat(outfile, cmd->operation);
    
    strcat(outfile, ".bw");

    fp = bwOpen(outfile, NULL, "w");
    bwCreateHdr(fp, cmd->binSize);
    fp->cl = bwCreateChromList(chrnames, chrlens, no_of_chrs);
    bwWriteHdr(fp);

    while (chr != NULL) {       
        if (chr->blacklist == 0 && chr->length > cmd->binSize) {
            if (chr->length > 10000000)
                printf("Writing: %s\n", chr->name);

            bwAddIntervalSpanSteps(fp, chr->name, start, cmd->binSize, cmd->binSize, &chr->coverages[curr->id][0], (uint32_t) 1);
            intervals = chr->coverages[curr->id] + 1;

            if(cmd->strand == -1 && strcmp(cmd->operation, INPUTS_ENDR) == 0) {
                
                for(i = 0; i <= chr->numberOfBins - 1; i++) {
                    if(intervals[i] > 0)
                        intervals[i] = -intervals[i]; 
                }
            }
            
            bwAppendIntervalSpanSteps(fp, intervals, (uint32_t) chr->numberOfBins - 1);
            
            if(cmd->strand == -1 && strcmp(cmd->operation, INPUTS_ENDR) == 0) {
                for(i = 0; i <= chr->numberOfBins - 1; i++) {
                    if(intervals[i] < 0) {
                        intervals[i] = -intervals[i];
                    }
                }
                
            }
        }

        chr = chr->next;
    }
    
    bwClose(fp);
    bwCleanup();

    if (chrnames) {
        for (i = 0; i < no_of_chrs; i++) {
            if (chrnames[i])
                free(chrnames[i]);
        }
        free(chrnames);
    }
    if (chrlens)
        free(chrlens);

    if (outfile)
        free(outfile);
}