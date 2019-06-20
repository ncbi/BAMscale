/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   scale.c
 * Author: pongorls
 * 
 * Created on November 28, 2018, 3:54 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libgen.h>

#include <bigWig.h>

#include "Definitions.h"
#include "main.h"
#include "scale.h"
#include "CHROMstruct.h"
#include "segmenter.h"
#include "binning.h"

float *scaleBins(float *carray, float scale, int nbins, float pseudocount) {
    int i = 0;
    
    for (i = 0; i < nbins; i++) {
        if (carray[i] > 0) {
            carray[i] = carray[i] * scale;
            carray[i] = roundf(carray[i] * 100.0) / 100.0;
        } else
            carray[i] = 0;
    }
    
    return carray;
}

void ScaleToSmallest(BAMFILES *head) {
    BAMFILES *curr = head;

    int smallestBAM = -1;

    while (curr != NULL) {
        if (curr->read_coverage != -1) {
            if (smallestBAM == -1) {
                smallestBAM = curr->read_coverage;
            }
            else {
                if (smallestBAM > curr->read_coverage && curr->read_coverage > 0) {
                    smallestBAM = curr->read_coverage;
                }
            }
        }

        curr = curr->next;
    }

    curr = head;

    while (curr != NULL) {
        if (curr->read_coverage > 0) {
            curr->scale = (float) curr->read_coverage / (float) smallestBAM;
            curr->scale = 1 / curr->scale;
        }

        curr = curr->next;
    }
}

void NoScale(BAMFILES *head) {
    BAMFILES *curr = head;

    while (curr != NULL) {
        if (curr->read_coverage != -1) {
            curr->scale = (float) 1;
        }

        curr = curr->next;
    }

    curr = head;
}

void ScaleToGenomeSize(BAMFILES *head, CHROMOSOMES *chead) {
    BAMFILES *curr = head;

    float genome_size = CalculateGenomeSize(chead) / 100;

    while (curr != NULL) {
        if (curr->read_coverage > 0) {
            curr->scale = 1 / ((float) curr->read_coverage / genome_size);
        }

        curr = curr->next;
    }
}

void ScaleGenomeCoverage(BAMFILES *head, CHROMOSOMES *chead) {
    float genome_size = CalculateGenomeSize(chead);
    BAMFILES *curr = head;

    while (curr != NULL) {
        curr->genome_scale = (float) 1 / (curr->base_coverage / genome_size);
        curr = curr->next;
    }
}

BAMFILES *ComputeSamplescales(BAMFILES *head, CHROMOSOMES *chead, int scale) {
    if (scale == 0) {
        NoScale(head);
    }

    if (scale == 1) {
        ScaleToSmallest(head);
    }

    if (scale == 2) {
        ScaleToGenomeSize(head, chead);
    }

    return head;
}

float *logTwoCoverageRatio(float *cov1, float *cov2, int nbins, float min_per_bin_cov) {
    int i = 0;
    float *carray = (float *) calloc(nbins + 1, sizeof (float));

    for (i = 0; i < nbins; i++) {
        if (cov1[i] >= min_per_bin_cov && cov2[i] >= min_per_bin_cov) {
            carray[i] = log2(((cov1[i]) / (cov2[i])));
            if (carray[i] > 1000)
                carray[i] = 1000;

            if (carray[i] < -1000)
                carray[i] = -1000;
        }
    }
       
    if (cov1)
        free(cov1);

    return carray;
}

float *OKseqRFD(float *cov1, float *cov2, int nbins, float min_per_bin_cov) {
    int i = 0;
    float *carray = (float *) calloc(nbins + 1, sizeof (float));

    for (i = 0; i < nbins; i++) {
        if (cov1[i] >= min_per_bin_cov || cov2[i] >= min_per_bin_cov) {
            carray[i] = ( (cov1[i]+0.0001) - (cov2[i]+0.0001)) / ((cov2[i]+0.0001) + (cov1[i]+0.0001));
        }
    }
       
    if (cov1)
        free(cov1);

    return carray;
}

float *SubtractCoverage(float *cov1, float *cov2, int nbins, float min_per_bin_cov) {
    int i = 0;
    float *carray = (float *) malloc(nbins * sizeof (float));

    for (i = 0; i < nbins; i++) {
        if (cov1[i] >= min_per_bin_cov || cov2[i] >= min_per_bin_cov) {
            carray[i] = cov1[i] - cov2[i];
        }
    }

    return carray;
}

float *CoverageRatio(float *cov1, float *cov2, int nbins, float min_per_bin_cov) {
    int i = 0;
    float *carray = (float *) malloc(nbins * sizeof (float));

    for (i = 0; i < nbins; i++) {
        if (cov1[i] >= min_per_bin_cov || cov2[i] >= min_per_bin_cov) {
            carray[i] = (cov1[i] / cov2[i]);
        }
    }

    return carray;
}

float *SignedCoverageRatio(float *cov1, float *cov2, int nbins, int min_per_bin_cov) {
    int i = 0;
    float *carray = (float *) malloc(nbins * sizeof (float));

    for (i = 0; i < nbins; i++) {
        if (cov1[i] >= min_per_bin_cov || cov2[i] >= min_per_bin_cov) {
            if (cov1[i] > cov2[i]) {
                carray[i] = (cov1[i] / cov2[i]);
            }
            else {
                carray[i] = -(cov2[i] / cov1[i]);
            }

        }
    }

    return carray;
}

CHRCOV *CalculateChromosomeRatio(CHROMOSOMES *curr, CHRCOV *chead, int s1, int s2, int ratioType, int min_per_bin_cov) {
    CHRCOV *ccurr = chead;

    if (ccurr == NULL) {
        chead = CreateChromCovStruct(curr->name, curr->id, curr->numberOfBins);
        ccurr = chead;
    }
    else {
        while (ccurr->next != NULL) {
            ccurr = ccurr->next;
        }

        ccurr->next = CreateChromCovStruct(curr->name, curr->id, curr->numberOfBins);
        ccurr = ccurr->next;
    }

    if (ratioType == 1) {
        ccurr->ratio = logTwoCoverageRatio(curr->coverages[s1], curr->coverages[s2], curr->numberOfBins, min_per_bin_cov);
    }
    else if (ratioType == 2) {
        ccurr->ratio = SubtractCoverage(curr->coverages[s1], curr->coverages[s2], curr->numberOfBins, min_per_bin_cov);
    }
    else if (ratioType == 3) {
        ccurr->ratio = CoverageRatio(curr->coverages[s1], curr->coverages[s2], curr->numberOfBins, min_per_bin_cov);
    }
    else if (ratioType == 4) {
        ccurr->ratio = SignedCoverageRatio(curr->coverages[s1], curr->coverages[s2], curr->numberOfBins, min_per_bin_cov);
    }
    
    return chead;
}

void PrintBedgraph(RATIOS *ptr, int binSize) {
    char outfile[250];
    int i = 0;
    CHRCOV *p = ptr->chrcovs;
    FILE * fp;

    strcpy(outfile, basename(ptr->sample1));
    strcat(outfile, "_vs_");
    strcat(outfile, basename(ptr->sample2));
    strcat(outfile, ".bedgraph");

    fp = fopen(outfile, "w+");

    while (p != NULL) {
        for (i = 0; i < p->nbins - 1; i++)
            fprintf(fp, "%s\t%d\t%d\t%.3f\n", p->name, i * binSize, (i + 1) * binSize, p->ratio[i]);
        p = p->next;
    }

    fclose(fp);
}

char *returnChrName(char *input) {
    char *ptr = strtok(input, "\t");

    while (ptr != NULL) {
        return ptr;
    }

    return NULL;
}

void PrintBedgraphOrdered(RATIOS *ptr, int binSize, char *chromfile) {
    FILE *handler = fopen(chromfile, "r");
    char line[BUFSIZ];
    char outfile[250];
    char *chrname = NULL;
    char *pos;
    int i = 0;
    CHRCOV *p = ptr->chrcovs;
    FILE * fp;

    strcpy(outfile, basename(ptr->sample1));
    strcat(outfile, "_vs_");
    strcat(outfile, basename(ptr->sample2));
    strcat(outfile, ".bedgraph");

    fp = fopen(outfile, "w+");

    while (fgets(line, sizeof (line), handler)) {
        if ((pos = strchr(line, '\n')) != NULL)
            *pos = '\0';

        chrname = returnChrName(line);
        p = ptr->chrcovs;

        if (chrname != NULL) {
            while (p != NULL) {
                if (strcmp(chrname, p->name) == 0) {
                    for (i = 0; i < p->nbins - 1; i++)
                        fprintf(fp, "%s\t%d\t%d\t%.3f\n", p->name, i * binSize, (i + 1) * binSize, p->ratio[i]);
                }

                p = p->next;
            }
        }
    }

    while (p != NULL) {
        for (i = 0; i < p->nbins - 1; i++)
            fprintf(fp, "%s\t%d\t%d\t%.3f\n", p->name, i * binSize, (i + 1) * binSize, p->ratio[i]);
        p = p->next;
    }

    fclose(fp);
    fclose(handler);
}

void PrintBigWigOrdered(RATIOS *ptr, int binSize, char *chromfile) {
    FILE *handler = NULL;
    char line[BUFSIZ];
    char outfile[250];
    char *pos, *chrname;
    int i = 0;
    char **chrnames = NULL;
    uint32_t *chrlens = NULL;
    uint32_t start = 0;
    uint32_t end = 0;
    int no_of_chrs = 0;
    CHRCOV *p = ptr->chrcovs;
    bigWigFile_t *fp;

    handler = fopen(chromfile, "r");

    while (fgets(line, sizeof (line), handler)) {
        if ((pos = strchr(line, '\n')) != NULL)
            *pos = '\0';

        chrname = returnChrName(line);
        p = ptr->chrcovs;

        if (chrname != NULL) {
            while (p != NULL) {
                if (strcmp(chrname, p->name) == 0) {
                    no_of_chrs++;
                }

                p = p->next;
            }
        }
    }

    fclose(handler);

    chrnames = (char **) malloc(no_of_chrs * sizeof (char *));
    chrlens = (uint32_t *) malloc(no_of_chrs * sizeof (uint32_t));
    no_of_chrs = 0;

    handler = fopen(chromfile, "r");

    while (fgets(line, sizeof (line), handler)) {
        if ((pos = strchr(line, '\n')) != NULL)
            *pos = '\0';

        chrname = returnChrName(line);
        p = ptr->chrcovs;

        if (chrname != NULL) {
            while (p != NULL) {
                if (strcmp(chrname, p->name) == 0) {
                    chrnames[no_of_chrs] = strdup(chrname);
                    chrlens[no_of_chrs] = (uint32_t) (p->nbins * binSize);
                    no_of_chrs++;
                }

                p = p->next;
            }
        }
    }

    fclose(handler);

    strcpy(outfile, basename(ptr->sample1));
    strcat(outfile, "_vs_");
    strcat(outfile, basename(ptr->sample2));
    strcat(outfile, ".bw");

    fp = bwOpen(outfile, NULL, "w");
    bwCreateHdr(fp, 10);
    fp->cl = bwCreateChromList(chrnames, chrlens, no_of_chrs);
    bwWriteHdr(fp);

    handler = fopen(chromfile, "r");

    while (fgets(line, sizeof (line), handler)) {
        if ((pos = strchr(line, '\n')) != NULL)
            *pos = '\0';

        chrname = returnChrName(line);
        p = ptr->chrcovs;

        if (chrname != NULL) {
            while (p != NULL) {
                if (strcmp(chrname, p->name) == 0) {
                    for (i = 0; i < p->nbins - 1; i++) {
                        start = (uint32_t) (i * binSize);
                        end = (uint32_t) ((i + 1) * binSize);
                        bwAddIntervals(fp, &p->name, &start, &end, &p->ratio[i], (uint32_t) 1);
                    }
                }

                p = p->next;
            }
        }
    }

    fclose(handler);
    bwClose(fp);
    bwCleanup();

    for (i = 0; i < no_of_chrs; i++) {
        if (chrnames[i]) {
            free(chrnames[i]);
        }
    }

    if (chrnames)
        free(chrnames);

    if (chrlens)
        free(chrlens);

}

RATIOS *CalculateRatiosAll(RATIOS *head, CHROMOSOMES *chead, BAMFILES *bhead, int no_of_samples, int min_per_bin_cov, int smoothbin, int binSize, char *chromsizes) {
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

        else
            PrintBedgraph(curr, binSize);

        bcurr = bcurr->next;
    }

    return head;
}