/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BEDstruct.c
 * Author: pongorls
 * 
 * Created on December 10, 2018, 8:02 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <pthread.h>
#include <htslib/sam.h>

#include "Definitions.h"
#include "main.h"
#include "BAMcoverage.h"
#include "multithreads.h"

char *BEDtoString(char *chr, int start, int end) {
    int char_written = 0;
    char *dest = NULL;

    if (start < 0)
        start = 0;

    dest = (char *) malloc(101 * sizeof (char));
    char_written = snprintf(dest, 100, "%s:%d-%d", chr, start, end);

    if (char_written >= 100 || char_written <= 3) {
        if (dest)
            free(dest);

        return NULL;
    }

    return dest;
}

int Read_filter_MultiCov(bam1_t *read, int paired_end) {
    if (paired_end == 1) {
        if (read->core.tid != read->core.mtid)
            return 0;

        if (abs((int) read->core.pos - (int) read->core.mpos) > 2000)
            return 0;
    }

    if ((int) read->core.qual < 10)
        return 0;

    if (read->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP))
        return 0;

    return 1;
}

void DeleteBEDs(PEAK *head) {
    PEAK *curr = head;

    while (head != NULL) {
        curr = head;
        head = head->next;

        if (curr->chr)
            free(curr->chr);

        if (curr->coord)
            free(curr->coord);

        if (curr->noSamples > 0) {
            if (curr->read_cov)
                free(curr->read_cov);

            if (curr->cov) {
                free(curr->cov);
            }

            if (curr->normalized)
                free(curr->normalized);
        }

        if (curr->end_str)
            free(curr->end_str);

        if (curr->start_str)
            free(curr->start_str);

        if (curr)
            free(curr);

    }
}

PEAK *CreateBEDentry() {
    PEAK *ptr = (PEAK *) malloc(sizeof (PEAK));

    ptr->next = NULL;
    ptr->binSize = -1;
    ptr->chr = NULL;
    ptr->cov = NULL;
    ptr->end = -1;
    ptr->end_str = NULL;
    ptr->id = -1;
    ptr->length = -1;
    ptr->coord = NULL;
    ptr->nbins = -1;
    ptr->start = -1;
    ptr->start_str = NULL;
    ptr->strand = 0;
    ptr->tid = -1;
    ptr->normalized = NULL;

    ptr->read_cov = NULL;
    ptr->noSamples = -1;

    return ptr;
}

void AllocateReadCovs(PEAK *head, int no_of_samples) {
    PEAK *curr = head;
    int i;

    while (curr != NULL) {
        curr->noSamples = no_of_samples;
        curr->read_cov = (int *) malloc(no_of_samples * sizeof (int));
        curr->normalized = (float *) malloc(no_of_samples * sizeof (float));

        for (i = 0; i < no_of_samples; i++) {
            curr->read_cov[i] = 0;
            curr->normalized[i] = 0.00;
        }

        curr = curr->next;
    }
}

void AllocateCovs(PEAK *head) {
    PEAK *curr = head;
    int i, j;

    while (curr != NULL) {
        if (curr->nbins > 0 && curr->noSamples > 0) {
            curr->cov = (float **) malloc(curr->noSamples * sizeof (float *));

            for (i = 0; i < curr->noSamples; i++) {
                curr->cov[i] = (float *) malloc(curr->nbins * sizeof (float));
                for (j = 0; j < curr->nbins; j++)
                    curr->cov[i][j] = 0.0;
            }
        }

        curr = curr->next;
    }
}

PEAK *AddBEDentry(PEAK *curr, char *BEDentry, int threadId) {
    char * ptr;
    int i = 0;

    if (curr == NULL) {
        curr = CreateBEDentry();
        curr->id = 0;
    } else {
        curr->next = CreateBEDentry();
        curr->next->id = curr->id + 1;
        curr = curr->next;
    }

    curr->tid = threadId;

    ptr = strtok(BEDentry, "\t");

    while (ptr != NULL) {
        if (i == 0) {
            curr->chr = strdup(ptr);
        }

        if (i == 1) {
            curr->start_str = strdup(ptr);
            curr->start = atoi(curr->start_str);
        }

        if (i == 2) {
            curr->end_str = strdup(ptr);
            curr->end = atoi(curr->end_str);
        }

        if (i == 5) {
            if (strncmp(ptr, "-", 1) == 0)
                curr->strand = -1;

            if (strncmp(ptr, "+", 1) == 0)
                curr->strand = 1;
        }

        i++;
        ptr = strtok(NULL, "\t");
    }

    if (curr->chr != NULL && curr->end_str != NULL && curr->start_str != NULL) {
        int len_str = strlen(curr->chr) + strlen(curr->start_str) + strlen(curr->end_str) + 2;
        curr->coord = (char *) malloc((len_str + 1) * sizeof (char));
        strcpy(curr->coord, curr->chr);
        strcat(curr->coord, ":");
        strcat(curr->coord, curr->start_str);
        strcat(curr->coord, "-");
        strcat(curr->coord, curr->end_str);

        curr->length = curr->end - curr->start;
    }

    return curr;
}

PEAK *ReadBED(char *BEDfilename, int nthreads) {
    FILE * handler = fopen(BEDfilename, "r");
    PEAK *head = NULL;
    PEAK *curr = head;
    char line[BUFSIZ];
    char *pos;
    int curr_thread = 0;

    while (fgets(line, sizeof (line), handler)) {
        if ((pos = strchr(line, '\n')) != NULL)
            *pos = '\0';

        if (head == NULL) {
            head = AddBEDentry(head, line, curr_thread);
            curr = head;
        } else
            curr = AddBEDentry(curr, line, curr_thread);

        curr_thread++;

        if (curr_thread >= nthreads)
            curr_thread = 0;
    }

    fclose(handler);
    return head;
}

void GetBEDCoveragesBAM(BAMFILES *bhead, PEAK *beds, CMDINPUT *cmd) {
    PEAK *bedhead = beds;
    BAMFILES *bcurr = bhead;

    while (bcurr != NULL) {
        samFile *fp_in = hts_open(bcurr->name, "r");
        hts_idx_t *idx = sam_index_load(fp_in, bcurr->name);
        bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
        bam1_t *aln = bam_init1();

        bedhead = beds;

        while (bedhead != NULL) {
            if (bedhead->coord) {
                hts_itr_t *iter = bam_itr_querys(idx, hdr, bedhead->coord);

                while (sam_itr_next(fp_in, iter, aln) > 0) {
                    if (Read_filter(aln, cmd) == 1) {
                        bedhead->read_cov[bcurr->id]++;
                    }
                }
                hts_itr_destroy(iter);
            }

            bedhead = bedhead->next;
        }

        bam_destroy1(aln);
        bam_hdr_destroy(hdr);
        hts_idx_destroy(idx);
        sam_close(fp_in);
        bcurr = bcurr->next;
    }
}

void *GetBEDFragmentCoveragesBAMmultithread(void * voidA) {
    THREADS *ptr = (THREADS *) voidA;
    char *filename = ptr->sample;
    PEAK *curr = ptr->phead;

    samFile *fp_in = hts_open(filename, "r"); //open bam file
    bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
    hts_idx_t *idx = sam_index_load(fp_in, filename);

    int read_strand = 0;
    int checkStrandInfo = 0;

    int frag_start = 0;
    int frag_end = 0;
    int aln_start = 0;
    int aln_end = 0;

    char *ext_coord = NULL;

    while (curr != NULL) {
        if (curr->tid == ptr->pid && curr->coord != NULL) {
            checkStrandInfo = 0;

            if (ptr->strand != 0) {
                if (curr->strand == 0) {
                    checkStrandInfo = 1;
                    printf("[ %d ] Error: %s has no strand info\n", ptr->pid, curr->coord);
                }
            }

            ext_coord = BEDtoString(curr->chr, curr->start - ptr->cmd->fragment_size, curr->end + ptr->cmd->fragment_size);

            if (checkStrandInfo == 0 && ext_coord != NULL) {
                hts_itr_t * iter = bam_itr_querys(idx, hdr, curr->coord);

                while (sam_itr_next(fp_in, iter, aln) > 0) {
                    if (Read_filter(aln, ptr->cmd) == 1) {
                        if (ptr->cmd->libtype == 1) {
                            aln_start = (int) aln->core.pos;
                            frag_start = (int) aln->core.mpos;
                            frag_end = (int) bam_endpos(aln);
                        } else {
                            aln_end = (int) bam_endpos(aln);

                            if (aln->core.flag & (BAM_FREVERSE)) {
                                frag_start = (int) aln->core.pos - ptr->cmd->fragment_size + ((int) bam_endpos(aln) - (int) aln->core.pos);
                                frag_end = (int) aln->core.pos;
                                aln_start = frag_end;
                            } else {
                                aln_start = (int) aln->core.pos + ptr->cmd->fragment_size;
                                frag_start = (int) aln->core.pos;
                                frag_end = aln_end + ptr->cmd->fragment_size;
                            }

                        }

                        if (aln_start > frag_start) {
                            if (frag_start <= curr->end && frag_end >= curr->start) {
                                if (ptr->strand != 0) {
                                    read_strand = ReadStrand(aln, ptr->paired_end);

                                    if (ptr->strand == 1) {
                                        if (curr->strand == read_strand) {
                                            curr->read_cov[ptr->sample_id]++;
                                        }
                                    }

                                    if (ptr->strand == -1) {
                                        if (curr->strand != read_strand) {
                                            curr->read_cov[ptr->sample_id]++;
                                        }
                                    }
                                } else
                                    curr->read_cov[ptr->sample_id]++;
                            }
                        }
                    }
                }
                hts_itr_destroy(iter);
            }

            if (ext_coord != NULL)
                free(ext_coord);
        }

        curr = curr->next;
    }

    bam_destroy1(aln);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(fp_in);

    return NULL;
}

void *GetBEDCoveragesBAMmultithread(void * voidA) {
    THREADS *ptr = (THREADS *) voidA;
    char *filename = ptr->sample;
    PEAK *curr = ptr->phead;

    samFile *fp_in = hts_open(filename, "r"); //open bam file
    bam_hdr_t *hdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
    hts_idx_t *idx = sam_index_load(fp_in, filename);

    int read_strand = 0;
    int checkStrandInfo = 0;

    while (curr != NULL) {
        if (curr->tid == ptr->pid && curr->coord != NULL) {
            checkStrandInfo = 0;

            if (ptr->strand != 0) {
                if (curr->strand == 0) {
                    checkStrandInfo = 1;
                    printf("[ %d ] Error: %s has no strand info\n", ptr->pid, curr->coord);
                }
            }

            if (checkStrandInfo == 0) {
                hts_itr_t *iter = bam_itr_querys(idx, hdr, curr->coord);

                while (sam_itr_next(fp_in, iter, aln) > 0) {
                    if (Read_filter(aln, ptr->cmd) == 1) {
                        if (ptr->strand != 0) {
                            read_strand = ReadStrand(aln, ptr->paired_end);

                            if (ptr->strand == 1) {
                                if (curr->strand == read_strand) {
                                    curr->read_cov[ptr->sample_id]++;
                                }
                            }

                            if (ptr->strand == -1) {
                                if (curr->strand != read_strand) {
                                    curr->read_cov[ptr->sample_id]++;
                                }
                            }
                        } else
                            curr->read_cov[ptr->sample_id]++;
                    }
                }
                hts_itr_destroy(iter);
            }
        }

        curr = curr->next;
    }

    bam_destroy1(aln);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(fp_in);

    return NULL;
}

void MultiCoverage(BAMFILES *bhead, PEAK *head, CMDINPUT *cmd) {
    BAMFILES *bcurr = bhead;
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
        threadStruct[i].chr = NULL;
        threadStruct[i].phead = NULL;
        threadStruct[i].cmd = cmd;
        threadStruct[i].bamfile = NULL;
        threadStruct[i].next = NULL;
    }

    while (bcurr != NULL) {
        printf("\t[ %d / %d ] %s\n", bcurr->id, cmd->no_of_samples - 1, bcurr->shortname);

        for (i = 0; i < cmd->threads; i++) {
            threadStruct[i].pid = i;
            if (threadStruct[i].chrname)
                free(threadStruct[i].chrname);
            threadStruct[i].sample = bcurr->name;
            threadStruct[i].paired_end = cmd->libtype;
            threadStruct[i].phead = head;
            threadStruct[i].sample_id = bcurr->id;
            threadStruct[i].strand = cmd->strand;
            threadStruct[i].cmd = cmd;
        }

        for (i = 0; i < cmd->threads; i++) {
            if (cmd->fragment_count_mode == 0)
                pthread_create(&thread_id[i], NULL, &GetBEDCoveragesBAMmultithread, (void *) &(threadStruct[i]));

            if (cmd->fragment_count_mode == 1)
                pthread_create(&thread_id[i], NULL, &GetBEDFragmentCoveragesBAMmultithread, (void *) &(threadStruct[i]));
        }

        for (i = 0; i < cmd->threads; i++)
            pthread_join(thread_id[i], NULL);

        bcurr = bcurr->next;
    }
    DestroyThreadStruct(&threadStruct, cmd->threads);
}

void CalculateFPKM(BAMFILES *bhead, PEAK *head) {
    BAMFILES *bcurr = bhead;
    PEAK *curr = head;

    float pm = 0;
    float kilo_length = 0;

    while (bcurr != NULL) {
        pm = (float) bcurr->read_coverage / 1000000;
        curr = head;

        while (curr != NULL) {
            curr->normalized[bcurr->id] = 0;
            
            if (curr->read_cov[bcurr->id] > 0) {
                kilo_length = (float) curr->length / 1000;
                curr->normalized[bcurr->id] = (float) curr->read_cov[bcurr->id] / (pm * kilo_length);
            }

            curr = curr->next;
        }

        bcurr = bcurr->next;
    }
}

void CalculateLibScaled(BAMFILES *bhead, PEAK *head) {
    BAMFILES *bcurr = bhead;
    PEAK *curr = head;

    while (bcurr != NULL) {
        curr = head;

        while (curr != NULL) {
            curr->normalized[bcurr->id] = 0;
            
            if (curr->read_cov[bcurr->id] > 0) {
                curr->normalized[bcurr->id] = (float) curr->read_cov[bcurr->id] * bcurr->scale;
            }

            curr = curr->next;
        }

        bcurr = bcurr->next;
    }
}

void CalculateTPM(BAMFILES *bhead, PEAK *head) {
    BAMFILES *bcurr = bhead;
    PEAK *curr = head;

    float kilo_length = 0;
    float sum_rpk = 0.0;

    while (bcurr != NULL) {
        curr = head;
        sum_rpk = 0.0;

        while (curr != NULL) {
            curr->normalized[bcurr->id] = 0;
            
            if (curr->read_cov[bcurr->id] > 0) {
                kilo_length = (float) curr->length / 1000;
                curr->normalized[bcurr->id] = (float) curr->read_cov[bcurr->id] / kilo_length;
                //printf("Sum %f\n", curr->normalized[bcurr->id]);
                sum_rpk += curr->normalized[bcurr->id];
            }

            curr = curr->next;
        }

        curr = head;
        sum_rpk = sum_rpk / 1000000;
        
        while (curr != NULL) {
            if (curr->read_cov[bcurr->id] > 0)
                curr->normalized[bcurr->id] = curr->normalized[bcurr->id] / sum_rpk;

            curr = curr->next;
        }

        bcurr = bcurr->next;
    }
}

void WriteMultiCovsRaw(BAMFILES *bhead, PEAK *head, int no_of_samples, char *outfile) {
    BAMFILES *bcurr = bhead;
    PEAK *curr = head;
    FILE *handler = fopen(outfile, "w");
    int i = 0;
    char * p = NULL;

    fprintf(handler, "coordinate");

    bcurr = bhead;

    while (bcurr != NULL) {
        p = strrchr(bcurr->name, '/');
        p = p ? p + 1 : (char *) bcurr->name;

        fprintf(handler, "\t%s", p);
        bcurr = bcurr->next;
    }

    fprintf(handler, "\n");

    while (curr != NULL) {
        fprintf(handler, "%s", curr->coord);

        for (i = 0; i < no_of_samples; i++) {
            fprintf(handler, "\t%d", curr->read_cov[i]);
        }

        fprintf(handler, "\n");
        curr = curr->next;
    }

    fclose(handler);

}

void WriteMultiCovsNormalized(BAMFILES *bhead, PEAK *head, int no_of_samples, char *outfile) {
    BAMFILES *bcurr = bhead;
    PEAK *curr = head;
    FILE *handler = fopen(outfile, "w");
    int i = 0;
    char * p = NULL;

    fprintf(handler, "coordinate");

    bcurr = bhead;

    while (bcurr != NULL) {
        p = strrchr(bcurr->name, '/');
        p = p ? p + 1 : (char *) bcurr->name;

        fprintf(handler, "\t%s", p);
        bcurr = bcurr->next;
    }

    fprintf(handler, "\n");

    while (curr != NULL) {
        fprintf(handler, "%s", curr->coord);

        for (i = 0; i < no_of_samples; i++) {
            fprintf(handler, "\t%f", curr->normalized[i]);
        }

        fprintf(handler, "\n");
        curr = curr->next;
    }

    fclose(handler);

}
