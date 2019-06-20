/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BAMcoverage.c
 * Author: pongorls
 * 
 * Created on November 28, 2018, 12:27 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <htslib/sam.h>
#include <pthread.h>
#include <bigWig.h>

#include "main.h"
#include "Definitions.h"
#include "BAMcoverage.h"
#include "CHROMstruct.h"
#include "BAMstructs.h"
#include "binning.h"
#include "scale.h"
#include "multithreads.h"
#include "Inputs.h"

int ReadStrand(bam1_t *read, int paired_end) {
    if (paired_end == 0) {
        if (read->core.flag & (BAM_FREVERSE)) {
            return 1;
        } else
            return -1;
    } else {
        if (read->core.flag & (BAM_FREAD1)) {
            if (read->core.flag & (BAM_FREVERSE)) {
                return 1;
            } else
                return -1;
        }

        if (read->core.flag & (BAM_FREAD2)) {
            if (read->core.flag & (BAM_FREVERSE)) {
                return -1;
            } else
                return 11;
        }
    }

    return 0;
}

int DetectLibraryType(BAMFILES *bhead) {
    samFile *fp_in = hts_open(bhead->name, "r");
    bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1();
    int counted = 0;
    int paired = 0;

    while (sam_read1(fp_in, hdr, aln) > 0 && counted < 1000000) {
        if (aln->core.mpos > -1)
            paired++;

        counted++;
    }

    bam_destroy1(aln);
    bam_hdr_destroy(hdr);
    sam_close(fp_in);

    if (paired > 0)
        return 1;

    return 0;
}

/*
 * 
 */
int Read_filter(bam1_t *read, CMDINPUT *cmd) { //int paired_end, int rmDuplicate, int rmProper, int rmUnmappedPair) {
    if (cmd->libtype == 1) {
        if (read->core.tid != read->core.mtid && cmd->filtDiffChr == 1)
            return 0;
        
        if(cmd->filtInsSize == 1) {
            if (abs((int) read->core.pos - (int) read->core.mpos) > cmd->max_insert_size || abs((int) read->core.pos - (int) read->core.mpos) < cmd->min_insert_size)
               return 0;
            }
        
        if (cmd->nounproper == 1 && !(read->core.flag & (BAM_FPROPER_PAIR)))
            return 0;
        
        if (cmd->remove_unmapped_pair == 1 && read->core.flag & (BAM_FUNMAP))
            return 0;
        
    }

    if ((int) read->core.qual < cmd->mapq)
        return 0;

    if (cmd->removeduplicates == 1 && read->core.flag & (BAM_FDUP))
        return 0;

    if (read->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL))
        return 0;

    return 1;
}

CHROMOSOMES *AddIDXcoverage(char *name, int coverage, int id, CHROMOSOMES *head) {
    CHROMOSOMES *curr = head;

    while (curr != NULL) {

        if (strcmp(name, curr->name) == 0) {
            curr->idxreads[id] += coverage;
            curr = NULL;
        }

        if (curr)
            curr = curr->next;
    }

    return head;
}

/**
 * Calculates reads from the index of the BAM file.
 * @param bamfile is the name of input BAM file (indexed)
 * @return <b>INT</b> showing the number of reads
 */
void GetChromosomeCoveragesIDX(CHROMOSOMES *head, BAMFILES *bhead) {
    int i;

    while (bhead != NULL) {
        samFile *fp_in = hts_open(bhead->name, "r");
        hts_idx_t *idx = sam_index_load(fp_in, bhead->name);
        bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
        for (i = 0; i < hdr->n_targets; i++) {
            uint64_t u, v;
            hts_idx_get_stat(idx, i, &u, &v);
            head = AddIDXcoverage(hdr->target_name[i], (int) u, bhead->id, head);
        }
        bam_hdr_destroy(hdr);
        hts_idx_destroy(idx);
        sam_close(fp_in);
        bhead = bhead->next;

    }
}

void GetGenomeCoveragesIDX(CHROMOSOMES *head, BAMFILES *bhead) {
    BAMFILES *bcurr = bhead;
    CHROMOSOMES *curr = head;

    while (bcurr != NULL) {
        bcurr->read_coverage = 0;
        bcurr = bcurr->next;
    }

    bcurr = bhead;

    while (bcurr != NULL) {
        curr = head;

        while (curr != NULL) {
            if (curr->idxreads && curr->blacklist == 0) {
                bcurr->read_coverage += curr->idxreads[bcurr->id];
            }

            curr = curr->next;
        }

        bcurr = bcurr->next;
    }
}

void CalculateCoverageOfReads(samFile *fp_in, hts_itr_t *iter, bam1_t *aln, int chrsize, char *chrname, CMDINPUT *cmd, BAMFILES *bamcurr) {
    while (sam_itr_next(fp_in, iter, aln) >= 0) {
        int tstrand = 0;
        
        if (Read_filter(aln, cmd) == 1) {
            tstrand = ReadStrand(aln, cmd->libtype);
            
            if(cmd->strand == 0 || tstrand == cmd->strand) {    
                if (cmd->fragment_count_mode == 1) {
                    if (cmd->libtype == 0)
                        bamcurr->read_coverage++;
                                   
                    else {
                        if (aln->core.flag & (BAM_FREAD1))
                            bamcurr->read_coverage++;
                    }
                } else {
                    bamcurr->read_coverage++;
                }    
            }
            
            else
                bamcurr->filtered_reads++;
            
        } else {
            bamcurr->filtered_reads++;
        }
    }
}

void *GetGenomeReadCoveragemultithread(void * voidA) {
    THREADS *ptr = (THREADS *) voidA;
    char *filename = ptr->sample;
    CHROMOSOMES *curr = ptr->chr;

    samFile *fp_in = hts_open(filename, "r"); //open bam file
    bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
    hts_idx_t *idx = sam_index_load(fp_in, filename);

    while (curr != NULL) {
        if (curr->tid == ptr->pid && curr->blacklist == 0) {
            if (curr->length > 10000000)
                fprintf(stderr, "\t\tThread [ %d ] is processing: %s\n", ptr->pid + 1, curr->name);
            hts_itr_t *iter = bam_itr_querys(idx, hdr, curr->name);
            CalculateCoverageOfReads(fp_in, iter, aln, curr->length, curr->name, ptr->cmd, ptr->bamfile);    
            hts_itr_destroy(iter);
            
        }
        curr = curr->next;
    }

    bam_destroy1(aln);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(fp_in);

    return NULL;
}

void MultiGenomeReadCoverage(CMDINPUT *cmd, CHROMOSOMES *chr) {
    BAMFILES *bamcurr = cmd->bamfiles;
    pthread_t thread_id[cmd->threads];
    int i = 0;

    fprintf(stderr, "\nImporting sequencing coverages\n");
    THREADS *threadStruct = (THREADS *) malloc(cmd->threads * sizeof (THREADS));
    for (i = 0; i < cmd->threads; i++) {
        threadStruct[i].pid = -1;
        threadStruct[i].chrname = NULL;
        threadStruct[i].sample = NULL;
        threadStruct[i].sample_id = -1;
        threadStruct[i].paired_end = 0;
        threadStruct[i].scale = 1.00;
        threadStruct[i].binSize = 0;
        threadStruct[i].pseudocount = 1;
        threadStruct[i].strand = 0;
        threadStruct[i].chr = chr;
        threadStruct[i].phead = NULL;
        threadStruct[i].cmd = cmd;
        threadStruct[i].bamfile = NULL;
        threadStruct[i].next = NULL;
    }

    while (bamcurr != NULL) {
        printf("\t[ %d / %d ] %s\n", bamcurr->id, cmd->no_of_samples - 1, bamcurr->shortname);

        for (i = 0; i < cmd->threads; i++) {
            threadStruct[i].pid = i;
            if (threadStruct[i].chrname)
                free(threadStruct[i].chrname);
            threadStruct[i].chrname = NULL;
            threadStruct[i].sample = bamcurr->name;
            threadStruct[i].sample_id = bamcurr->id;
            threadStruct[i].bamfile = bamcurr;
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_create(&thread_id[i], NULL, &GetGenomeReadCoveragemultithread, (void *) &(threadStruct[i]));
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_join(thread_id[i], NULL);
        }

        bamcurr = bamcurr->next;
    }
    DestroyThreadStruct(&threadStruct, cmd->threads);
}

/**
 * Calculates reads from the index of the BAM file.
 * @param bamfile is the name of input BAM file (indexed)
 * @return <b>INT</b> showing the number of reads
 */
void GetChromosomeCoveragesBAM(CHROMOSOMES *head, BAMFILES *bhead, CMDINPUT *cmd) {
    CHROMOSOMES *curr = head;

    while (bhead != NULL) {
        samFile *fp_in = hts_open(bhead->name, "r");
        hts_idx_t *idx = sam_index_load(fp_in, bhead->name);
        bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
        hts_itr_t *iter = NULL;
        bam1_t *aln = bam_init1();

        curr = head;

        while (curr != NULL) {
            if (curr->blacklist == 0) {
                //printf("Getting chromosome coverage for chr %s in sample %d\n", curr->name, bhead->id);
                curr->idxreads[bhead->id] = 0;
                iter = bam_itr_querys(idx, hdr, curr->name);

                while (sam_itr_next(fp_in, iter, aln) > 0) {
                    if (Read_filter(aln, cmd) == 1) {
                        curr->idxreads[bhead->id]++;
                    }
                }
            }

            curr = curr->next;
        }

        bam_destroy1(aln);
        bam_hdr_destroy(hdr);
        hts_idx_destroy(idx);
        sam_close(fp_in);
        bhead = bhead->next;
    }
}

char *BEDentryToCoord(char *input) {
    char *coord = NULL;
    char *ncoord = NULL;
    char *ptr;
    int i = 0;

    ptr = strtok(input, "\t");

    while (ptr != NULL) {
        if (i == 0) {
            ncoord = (char *) malloc(sizeof (char) * (strlen(ptr) + 2));
            strcpy(ncoord, ptr);
            strcat(ncoord, ":");
        }

        if (i == 1) {
            ncoord = (char *) realloc(coord, (strlen(coord) + strlen(ptr) + 2) * sizeof (char));
            strcat(ncoord, ptr);
            strcat(ncoord, "-");
        }

        if (i == 2) {
            ncoord = (char *) realloc(coord, (strlen(coord) + strlen(ptr) + 1) * sizeof (char));
            strcat(ncoord, ptr);
        }

        if (ncoord != NULL)
            coord = ncoord;

        ncoord = NULL;
        i++;
        ptr = strtok(NULL, "\t");
    }

    return coord;
}

char *BEDentryChr(char *input) {
    char *coord = NULL;
    char *ptr;
    int i = 0;

    ptr = strtok(input, "\t");

    while (ptr != NULL) {
        if (i == 0) {
            coord = strdup(ptr);
        }

        i++;
        ptr = strtok(NULL, "\t");
    }

    return coord;
}

void SubtractBlacklistedBEDS(char *filename, CHROMOSOMES *head, BAMFILES *bhead, int paired_end) {
    char line[BUFSIZ];
    char linecpy[BUFSIZ];
    char *coord, *chrname;
    int count = 0;

    while (bhead != NULL) {
        samFile *fp_in = hts_open(bhead->name, "r");
        hts_idx_t *idx = sam_index_load(fp_in, bhead->name);
        bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
        hts_itr_t *iter = NULL;
        bam1_t *aln = bam_init1();
        FILE * handler = fopen(filename, "r");

        while (fgets(line, sizeof (line), handler)) {
            if (line[strlen(line) - 1] == '\n')
                line[strlen(line) - 1] = '\0';

            strcpy(linecpy, line);

            coord = BEDentryToCoord(line);
            chrname = BEDentryChr(linecpy);

            count = 0;

            iter = bam_itr_querys(idx, hdr, coord);

            while (sam_itr_next(fp_in, iter, aln) > 0) {
                bhead->read_coverage--;
                count--;
            }

            head = AddIDXcoverage(chrname, count, bhead->id, head);


            if (coord)
                free(coord);

            if (chrname)
                free(chrname);
        }

        bam_destroy1(aln);
        bam_hdr_destroy(hdr);
        hts_idx_destroy(idx);
        sam_close(fp_in);
        fclose(handler);

        bhead = bhead->next;
    }
}

int *CalculateCoverage(samFile *fp_in, hts_itr_t *iter, bam1_t *aln, int chrsize, char *chrname, CMDINPUT *cmd, BAMFILES *bamcurr) {
    int *chr_cov = (int *) calloc(chrsize+1, sizeof (int));
    int i = 0;
    int j = 0;
    
    if (chr_cov == NULL) {
        printf("ERROR: could not allocate memory for coverage track of %s\n", chrname);
        FreeAllocatedData();
        exit(0);
    }

    while (sam_itr_next(fp_in, iter, aln) >= 0) {
        int tstrand = 0;
        int pos = -1; //left most position of alignment in zero based coordianate (+1)
        int end = -1; //end of read/fragment
        int ins_start = -1;
        int ins_end = -1;
        
        if (Read_filter(aln, cmd) == 1) {
            tstrand = ReadStrand(aln, cmd->libtype);
            
            if(cmd->strand == 0 || tstrand == cmd->strand) {    
                //printf("next\n");
                uint32_t *cigar = bam_get_cigar(aln);
                int fragend = chrsize;
                pos = (int) aln->core.pos;
                end = (int) bam_endpos(aln);

                //check straqndedness??
                if(cmd->strand != 0) {
                    if(tstrand != cmd->strand) {
                        pos = -1;
                    }
                }

                if (cmd->fragment_count_mode == 1) {
                    if (cmd->libtype == 0) {
                        bamcurr->read_coverage++;
                        if (tstrand == 1) {
                            ins_start = end;
                            ins_end = ins_start + cmd->fragment_size - (end - pos);
                        } else {
                            ins_start = end - cmd->fragment_size;
                            ins_end = pos;
                        }                   
                    } else {
                        if (aln->core.flag & (BAM_FREAD1))
                            bamcurr->read_coverage++;

                        if ((int) aln->core.pos < (int) aln->core.mpos) {
                            ins_start = end;
                            ins_end = (int) aln->core.mpos;
                        } else if ((int) aln->core.pos > (int) aln->core.mpos) {
                            ins_start = -1;
                            ins_end = -1;
                        } else {
                            if (aln->core.flag & (BAM_FREAD1)) {
                                ins_start = -1;
                                ins_end = -1;
                            }
                        }
                    }
                } else {
                    ins_start = -1;
                    ins_end = -1;
                    bamcurr->read_coverage++;
                }

                for(i = 0; i < aln->core.n_cigar; i++) {
                    int cop = cigar[i] & BAM_CIGAR_MASK;
                    int cl = cigar[i] >> BAM_CIGAR_SHIFT;
                    int endpos = pos;
                    
                    switch (cop) {
                        case BAM_CMATCH:
                            endpos = pos + cl;
                            break;

                        case BAM_CHARD_CLIP:
                            break;

                        case BAM_CSOFT_CLIP:
                            break;

                        case BAM_CDEL: // deletion, no coverage
                            pos+=cl;
                            break;

                        case BAM_CPAD:
                            pos+=cl;
                            break;

                        case BAM_CINS:
                            break;

                        case BAM_CREF_SKIP: // splicing bases, , no coverage
                            pos+=cl;
                            break;

                        default:
                            fprintf(stderr, "Unhandled cigar_op %d:%d\n", cop, cl);
                            printf("?");
                    }
                    
                    if(endpos > pos) {
                        for(j = pos; j < endpos; j++) {
                            if(j < fragend && j < chrsize && j >= 0)
                                chr_cov[j]++;
                        }
                        
                        bamcurr->base_coverage += (endpos - pos);
                        
                        pos = endpos;
                    }
                }

                if (ins_start != -1 && ins_end != -1) {
                    if (ins_start < 0)
                        ins_start = 0;

                    if (ins_end >= chrsize)
                        ins_end = chrsize - 1;

                    bamcurr->base_coverage += (ins_start - ins_end);

                    for (i = ins_start; i <= ins_end; i++) {
                        chr_cov[i]++;
                    }
                }
            } else {
                    bamcurr->read_coverage++;
                }
        } else {
            bamcurr->filtered_reads++;
        }
    }
    
    return chr_cov;
}


void GetGenomeCoverageRNA(CMDINPUT *cmd, CHROMOSOMES *head, char *outfile) {
    char *filename = cmd->bamfiles->name;

    samFile *fp_in = hts_open(filename, "r"); //open bam file
    bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
    hts_idx_t *idx = sam_index_load(fp_in, filename);
    CHROMOSOMES *curr = head;
    bigWigFile_t *fp = NULL;
    int no_of_chrs = CountNumberOfChromosomes(cmd->chr);
    float *average = (float *)calloc(1, sizeof(float));
    
    char **chrnames = GetChromosomeNames(cmd->chr, no_of_chrs);
    uint32_t *chrlens = GetChrLens(cmd->chr, no_of_chrs);
    int invert = 0;
    uint32_t start = 0;
    
    if(cmd->strand == -1 && strcmp(cmd->operation, INPUTS_RSTRRNA) == 0)
        invert = 1;
    
    fp = bwOpen(outfile, NULL, "w");
    bwCreateHdr(fp, 7);
    fp->cl = bwCreateChromList(chrnames, chrlens, no_of_chrs);
    bwWriteHdr(fp);
    
    while (curr != NULL) {
        if (curr->blacklist == 0) {
            if (curr->length > 10000000)
                fprintf(stderr, "\t\tProcessing: %s\n",curr->name);

            hts_itr_t *iter = bam_itr_querys(idx, hdr, curr->name);
            
            int *coverage = CalculateCoverage(fp_in, iter, aln, curr->length, curr->name, cmd, cmd->bamfiles);
            int *binres = (int *) calloc(curr->length+1, sizeof (int));
            int currbin = 0;
            int binsize = cmd->binSize;
            int maxdiff = 5;
                        
            for(int i = 0; i < curr->length - 1; i++) {
                currbin = i / binsize;
                
                if(abs(coverage[i] - coverage[i+1]) > maxdiff) {
                    //if(coverage[i] < 10 || coverage[i+1] < 10) {
                        if(currbin == (i+1) / binsize) {
                            int end = currbin*cmd->binSize+cmd->binSize;
                            
                            if(end > curr->length)
                                end = curr->length;
                            
                            for(int j = currbin*cmd->binSize; j < end; j++) {
                                //printf("%d\t%d\n", j, curr->length);
                                binres[j] = 1;
                            }
                        }
                    //}
                }
            }
                       
            int prev_type = 0;
            start = 0;
            
            for(int i = 0; i < curr->numberOfBins-1; i++) {
                start = (uint32_t) (i * cmd->binSize);
                
                if(binres[i*cmd->binSize] == 0) {
                    average[0] = 0;
                    
                    int end = (i+1) * cmd->binSize;
                            
                    if(end > curr->length)
                        end = curr->length;
                    
                    for(int j = i*cmd->binSize; j < end; j++)
                        average[0] += (float)coverage[j];

                    if(average[0] > 0)
                        average[0] = average[0] / (float) cmd->binSize;

                    average[0] = average[0] * cmd->bamfiles->genome_scale;
                    
                    if(cmd->strand == -1 && invert == 1 && average[0] > 0)
                        average[0] = -average[0];

                    if(i == 0) {
                        bwAddIntervalSpanSteps(fp, curr->name, start, cmd->binSize, cmd->binSize, &average[0], (uint32_t) 1);
                    }
                    
                    else {
                        if(prev_type != 0) {
                            start = (uint32_t) (i*cmd->binSize);
                            bwAddIntervalSpanSteps(fp, curr->name, start, cmd->binSize, cmd->binSize, &average[0], (uint32_t) 1);
                        }
                        
                        else {
                            bwAppendIntervalSpanSteps(fp, &average[0], (uint32_t) 1);
                        }
                    }
                    
                    prev_type = 0;
                }
                
                else {
                    int j = 0;
                    
                    if(prev_type == 0) {
                        j = 1;
                        
                        average[0] = (float)coverage[i*cmd->binSize];
                        average[0] = average[0] * cmd->bamfiles->genome_scale;
                        
                        if(cmd->strand == -1 && invert == 1 && average[0] > 0)
                            average[0] = -average[0];
                        
                        bwAddIntervalSpanSteps(fp, curr->name, start, 1, 1, &average[0], (uint32_t) 1);
                    }
                    
                    int end = (i + 1)*cmd->binSize;
                            
                    if(end > curr->length)
                        end = curr->length;
                    
                    for(j = j + i*cmd->binSize; j < end; j++) {
                        average[0] = (float)coverage[j];
                        average[0] = average[0] * cmd->bamfiles->genome_scale;
                        
                        if(cmd->strand == -1 && invert == 1 && average[0] > 0)
                            average[0] = -average[0];
                        
                        bwAppendIntervalSpanSteps(fp, &average[0], (uint32_t) 1);
                    }
                    
                    prev_type = 1;
                }
            }
            
            hts_itr_destroy(iter);
            free(coverage);
            free(binres);
        }
        curr = curr->next;
    }
    
    bwClose(fp);
    bwCleanup();
    bam_destroy1(aln);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(fp_in);
    free(chrlens);
    free(average);
    
    for(int i = 0; i <= no_of_chrs; i++) {
        if(chrnames[i])
            free(chrnames[i]);
    }

    if(chrnames)
        free(chrnames);

    if(outfile)
        free(outfile);
}

void *GetGenomeCoveragemultithread(void * voidA) {
    THREADS *ptr = (THREADS *) voidA;
    char *filename = ptr->sample;
    CHROMOSOMES *curr = ptr->chr;

    samFile *fp_in = hts_open(filename, "r"); //open bam file
    bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
    hts_idx_t *idx = sam_index_load(fp_in, filename);

    while (curr != NULL) {
        if (curr->tid == ptr->pid && curr->blacklist == 0) {
            if (curr->length > 10000000)
                fprintf(stderr, "\t\tThread [ %d ] is processing: %s\n", ptr->pid + 1, curr->name);

            hts_itr_t *iter = bam_itr_querys(idx, hdr, curr->name);
            
            int *coverage = CalculateCoverage(fp_in, iter, aln, curr->length, curr->name, ptr->cmd, ptr->bamfile);

            hts_itr_destroy(iter);
            free(coverage);

        }
        curr = curr->next;
    }

    bam_destroy1(aln);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(fp_in);

    return NULL;
}

void MultiGenomeCoverage(CMDINPUT *cmd, CHROMOSOMES *chr) {
    BAMFILES *bamcurr = cmd->bamfiles;
    pthread_t thread_id[cmd->threads];
    int i = 0;

    fprintf(stderr, "\nImporting sequencing coverages\n");
    THREADS *threadStruct = (THREADS *) malloc(cmd->threads * sizeof (THREADS));
    for (i = 0; i < cmd->threads; i++) {
        threadStruct[i].pid = -1;
        threadStruct[i].chrname = NULL;
        threadStruct[i].sample = NULL;
        threadStruct[i].sample_id = -1;
        threadStruct[i].paired_end = 0;
        threadStruct[i].scale = 1.00;
        threadStruct[i].binSize = 0;
        threadStruct[i].pseudocount = 1;
        threadStruct[i].strand = 0;
        threadStruct[i].chr = chr;
        threadStruct[i].phead = NULL;
        threadStruct[i].cmd = cmd;
        threadStruct[i].bamfile = NULL;
        threadStruct[i].next = NULL;
    }

    while (bamcurr != NULL) {
        printf("\t[ %d / %d ] %s\n", bamcurr->id, cmd->no_of_samples - 1, bamcurr->shortname);

        for (i = 0; i < cmd->threads; i++) {
            threadStruct[i].pid = i;
            if (threadStruct[i].chrname)
                free(threadStruct[i].chrname);
            
            threadStruct[i].chrname = NULL;
            threadStruct[i].sample = bamcurr->name;
            threadStruct[i].sample_id = bamcurr->id;
            threadStruct[i].bamfile = bamcurr;
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_create(&thread_id[i], NULL, &GetGenomeCoveragemultithread, (void *) &(threadStruct[i]));
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_join(thread_id[i], NULL);
        }

        bamcurr = bamcurr->next;
        
        if(cmd->strandsplit)
            cmd->strand = -1;
    }
    
    DestroyThreadStruct(&threadStruct, cmd->threads);
}

void *GetGenomeBaseCoveragemultithread(void * voidA) {
    THREADS *ptr = (THREADS *) voidA;
    char *filename = ptr->sample;
    CHROMOSOMES *curr = ptr->chr;

    samFile *fp_in = hts_open(filename, "r"); //open bam file
    bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
    hts_idx_t *idx = sam_index_load(fp_in, filename);

    while (curr != NULL) {
        if (curr->tid == ptr->pid && curr->blacklist == 0) {
            if (curr->length > 10000000)
                fprintf(stderr, "\t\tThread [ %d ] is processing: %s\n", ptr->pid + 1, curr->name);

            hts_itr_t *iter = bam_itr_querys(idx, hdr, curr->name);
            
            int *coverage = CalculateCoverage(fp_in, iter, aln, curr->length, curr->name, ptr->cmd, ptr->bamfile);

            if (curr->coverages[ptr->sample_id])
                free(curr->coverages[ptr->sample_id]);
            
            curr->coverages[ptr->sample_id] = BinCoverage(coverage, curr->length, ptr->cmd->binSize, curr->numberOfBins);

            hts_itr_destroy(iter);
            free(coverage);

        }
        curr = curr->next;
    }

    bam_destroy1(aln);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(fp_in);

    return NULL;
}

void MultiGenomeBaseCoverage(CMDINPUT *cmd, CHROMOSOMES *chr) {
    BAMFILES *bamcurr = cmd->bamfiles;
    pthread_t thread_id[cmd->threads];
    int i = 0;

    fprintf(stderr, "\nImporting sequencing coverages\n");
    THREADS *threadStruct = (THREADS *) malloc(cmd->threads * sizeof (THREADS));
    for (i = 0; i < cmd->threads; i++) {
        threadStruct[i].pid = -1;
        threadStruct[i].chrname = NULL;
        threadStruct[i].sample = NULL;
        threadStruct[i].sample_id = -1;
        threadStruct[i].paired_end = 0;
        threadStruct[i].scale = 1.00;
        threadStruct[i].binSize = 0;
        threadStruct[i].pseudocount = 1;
        threadStruct[i].strand = 0;
        threadStruct[i].chr = chr;
        threadStruct[i].phead = NULL;
        threadStruct[i].cmd = cmd;
        threadStruct[i].bamfile = NULL;
        threadStruct[i].next = NULL;
    }

    while (bamcurr != NULL) {
        printf("\t[ %d / %d ] %s\n", bamcurr->id, cmd->no_of_samples - 1, bamcurr->shortname);

        for (i = 0; i < cmd->threads; i++) {
            threadStruct[i].pid = i;
            if (threadStruct[i].chrname)
                free(threadStruct[i].chrname);
            
            threadStruct[i].chrname = NULL;
            threadStruct[i].sample = bamcurr->name;
            threadStruct[i].sample_id = bamcurr->id;
            threadStruct[i].bamfile = bamcurr;
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_create(&thread_id[i], NULL, &GetGenomeBaseCoveragemultithread, (void *) &(threadStruct[i]));
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_join(thread_id[i], NULL);
        }

        bamcurr = bamcurr->next;
        
        if(cmd->strandsplit)
            cmd->strand = -1;
    }
    
    DestroyThreadStruct(&threadStruct, cmd->threads);
}

void *ScaleBinsmultithread(void * voidA) {
    THREADS *ptr = (THREADS *) voidA;
    CHROMOSOMES *curr = ptr->chr;

    while (curr != NULL) {
        if (curr->tid == ptr->pid && curr->blacklist == 0) {
            curr->coverages[ptr->sample_id] = scaleBins(curr->coverages[ptr->sample_id], ptr->bamfile->genome_scale, curr->numberOfBins, (float) 1);
        }

        curr = curr->next;
    }

    return NULL;
}

void MultiGenomeScaler(CMDINPUT *cmd, CHROMOSOMES *chr) {
    BAMFILES *bamcurr = cmd->bamfiles;
    pthread_t thread_id[cmd->threads];
    int i = 0;

    THREADS *threadStruct = (THREADS *) malloc(cmd->threads * sizeof (THREADS));
    for (i = 0; i < cmd->threads; i++) {
        threadStruct[i].pid = -1;
        threadStruct[i].chrname = NULL;
        threadStruct[i].sample = NULL;
        threadStruct[i].sample_id = -1;
        threadStruct[i].paired_end = 0;
        threadStruct[i].scale = 1.00;
        threadStruct[i].binSize = 0;
        threadStruct[i].pseudocount = 1;
        threadStruct[i].strand = 0;
        threadStruct[i].chr = chr;
        threadStruct[i].phead = NULL;
        threadStruct[i].cmd = cmd;
        threadStruct[i].bamfile = NULL;
        threadStruct[i].next = NULL;
    }

    while (bamcurr != NULL) {
        printf("\t[ %d / %d ] %s\n", bamcurr->id, cmd->no_of_samples - 1, bamcurr->shortname);

        for (i = 0; i < cmd->threads; i++) {
            threadStruct[i].pid = i;
            threadStruct[i].sample = bamcurr->name;
            threadStruct[i].sample_id = bamcurr->id;
            threadStruct[i].bamfile = bamcurr;
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_create(&thread_id[i], NULL, &ScaleBinsmultithread, (void *) &(threadStruct[i]));
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_join(thread_id[i], NULL);
        }

        bamcurr = bamcurr->next;
    }
    DestroyThreadStruct(&threadStruct, cmd->threads);
}

void *SmoothBinsmultithread(void * voidA) {
    THREADS *ptr = (THREADS *) voidA;
    CHROMOSOMES *curr = ptr->chr;

    while (curr != NULL) {
        if (curr->tid == ptr->pid && curr->blacklist == 0) {
            curr->coverages[ptr->sample_id] = QuicksmoothenBins(curr->coverages[ptr->sample_id], ptr->cmd->smoothBin, curr->numberOfBins); //was 200
        }

        curr = curr->next;
    }
    return NULL;
}

void MultiGenomeSmoother(CMDINPUT *cmd, CHROMOSOMES *chr) {
    BAMFILES *bamcurr = cmd->bamfiles;
    pthread_t thread_id[cmd->threads];
    int i = 0;

    THREADS *threadStruct = (THREADS *) malloc(cmd->threads * sizeof (THREADS));
    for (i = 0; i < cmd->threads; i++) {
        threadStruct[i].pid = -1;
        threadStruct[i].chrname = NULL;
        threadStruct[i].sample = NULL;
        threadStruct[i].sample_id = -1;
        threadStruct[i].paired_end = 0;
        threadStruct[i].scale = 1.00;
        threadStruct[i].binSize = 0;
        threadStruct[i].pseudocount = 1;
        threadStruct[i].strand = 0;
        threadStruct[i].chr = chr;
        threadStruct[i].phead = NULL;
        threadStruct[i].cmd = cmd;
        threadStruct[i].bamfile = NULL;
        threadStruct[i].next = NULL;
    }

    while (bamcurr != NULL) {
        printf("\t[ %d / %d ] %s\n", bamcurr->id, cmd->no_of_samples - 1, bamcurr->shortname);

        for (i = 0; i < cmd->threads; i++) {
            threadStruct[i].pid = i;
            threadStruct[i].sample = bamcurr->name;
            threadStruct[i].sample_id = bamcurr->id;
            threadStruct[i].bamfile = bamcurr;
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_create(&thread_id[i], NULL, &SmoothBinsmultithread, (void *) &(threadStruct[i]));
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_join(thread_id[i], NULL);
        }

        bamcurr = bamcurr->next;
    }
    DestroyThreadStruct(&threadStruct, cmd->threads);
}

void *TransformBinsmultithread(void * voidA) {
    THREADS *ptr = (THREADS *) voidA;
    CHROMOSOMES *curr = ptr->chr;

    while (curr != NULL) {
        if (curr->tid == ptr->pid && curr->blacklist == 0) {
            if (strcmp(ptr->cmd->operation, INPUTS_LOG2) == 0 || strcmp(ptr->cmd->operation, INPUTS_ENDR) == 0 || strcmp(ptr->cmd->operation, INPUTS_END) == 0 || strcmp(ptr->cmd->operation, INPUTS_REP) == 0)
                curr->coverages[ptr->sample_id] = logTwoCoverageRatio(curr->coverages[ptr->sample_id], curr->coverages[0], curr->numberOfBins, 0.1);

            if (strcmp(ptr->cmd->operation, INPUTS_RATIO) == 0)
                curr->coverages[ptr->sample_id] = CoverageRatio(curr->coverages[ptr->sample_id], curr->coverages[0], curr->numberOfBins, 0.1);

            if (strcmp(ptr->cmd->operation, INPUTS_SUBSTRACT) == 0)
                curr->coverages[ptr->sample_id] = SubtractCoverage(curr->coverages[ptr->sample_id], curr->coverages[0], curr->numberOfBins, 0.1);
            
            if (strcmp(ptr->cmd->operation, INPUTS_RFD) == 0)
                curr->coverages[ptr->sample_id] = OKseqRFD(curr->coverages[ptr->sample_id], curr->coverages[0], curr->numberOfBins, 0.1);
        }

        curr = curr->next;
    }

    return NULL;
}

void MultiGenomeTransform(CMDINPUT *cmd, CHROMOSOMES *chr) {
    BAMFILES *bamcurr = cmd->bamfiles->next;
    pthread_t thread_id[cmd->threads];
    int i = 0;

    THREADS *threadStruct = (THREADS *) malloc(cmd->threads * sizeof (THREADS));
    for (i = 0; i < cmd->threads; i++) {
        threadStruct[i].pid = -1;
        threadStruct[i].chrname = NULL;
        threadStruct[i].sample = NULL;
        threadStruct[i].sample_id = -1;
        threadStruct[i].paired_end = 0;
        threadStruct[i].scale = 1.00;
        threadStruct[i].binSize = 0;
        threadStruct[i].pseudocount = 1;
        threadStruct[i].strand = 0;
        threadStruct[i].chr = chr;
        threadStruct[i].phead = NULL;
        threadStruct[i].cmd = cmd;
        threadStruct[i].bamfile = NULL;
        threadStruct[i].next = NULL;
    }

    while (bamcurr != NULL) {
        printf("\t[ %d / %d ] %s\n", bamcurr->id, cmd->no_of_samples - 1, bamcurr->shortname);

        for (i = 0; i < cmd->threads; i++) {
            threadStruct[i].pid = i;
            threadStruct[i].sample = bamcurr->name;
            threadStruct[i].sample_id = bamcurr->id;
            threadStruct[i].bamfile = bamcurr;
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_create(&thread_id[i], NULL, &TransformBinsmultithread, (void *) &(threadStruct[i]));
        }

        for (i = 0; i < cmd->threads; i++) {
            pthread_join(thread_id[i], NULL);
        }

        bamcurr = bamcurr->next;
    }
    DestroyThreadStruct(&threadStruct, cmd->threads);
}

/*
 *
 */
void CalculateCoverageOfChromosomeBins(CHROMOSOMES *head, BAMFILES *bhead, int paired_end, int bin_size, int pseudocount, CMDINPUT *cmd) {
    BAMFILES *bcurr = bhead;
    CHROMOSOMES *curr = head;
    int *cov = NULL;

    while (bcurr != NULL) {
        curr = head;
        samFile *fp_in = hts_open(bcurr->name, "r"); //open bam file
        bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
        bam1_t *aln = bam_init1(); //initialize an alignment
        hts_idx_t *idx = sam_index_load(fp_in, bcurr->name);

        while (curr != NULL) {
            if (curr->blacklist == 0) {
                printf("\tChr: %s, sample [ %d ] and %d\n", curr->name, bcurr->id, curr->blacklist);
                hts_itr_t *iter = bam_itr_querys(idx, bamHdr, curr->name);
                curr->coverages[bcurr->id] = BinCoverage(cov, curr->length, bin_size, curr->numberOfBins);
                curr->coverages[bcurr->id] = scaleBins(curr->coverages[bcurr->id], bcurr->scale, curr->numberOfBins, (float) pseudocount);
                
                hts_itr_destroy(iter);
                if (cov)
                    free(cov);

                cov = NULL;
            }

            curr = curr->next;
        }

        hts_idx_destroy(idx);
        bam_destroy1(aln);
        sam_close(fp_in);

        bcurr = bcurr->next;
    }
}

void SmoothenAllChromosomeBins(CHROMOSOMES *head, BAMFILES *bhead, int smoothBinNum) {
    BAMFILES *bcurr = bhead;
    CHROMOSOMES *curr = head;

    printf("Calculating smoothened Bins\n");

    while (bcurr != NULL) {
        curr = head;

        while (curr != NULL) {
            if (curr->blacklist == 0) {
                printf("\tChr: %s, sample [ %d ]\n", curr->name, bcurr->id);
                curr->coverages[bcurr->id] = smoothenBins(&curr->coverages[bcurr->id], smoothBinNum, curr->numberOfBins);
            }

            curr = curr->next;
        }

        bcurr = bcurr->next;
    }
}
