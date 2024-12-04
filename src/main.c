/* 
 * File:   main.c
 * Author: pongorls
 *
 * Created on November 28, 2018, 11:55 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>
#include <ctype.h>
#include <getopt.h>
#include <inttypes.h>

#include <bigWig.h>

#include "main.h"
#include "Definitions.h"
#include "BAMstructs.h"
#include "CHROMstruct.h"
#include "BAMcoverage.h"
#include "scale.h"
#include "segmenter.h"
#include "multithreads.h"
#include "BEDstruct.h"
#include "Inputs.h"
#include "Writer.h"
#include <htslib/sam.h>

BAMFILES *BAMhead = NULL;
BAMFILES *BAMcurr = NULL;

CHROMOSOMES *CHROMhead = NULL;
CHROMOSOMES *CHROMcurr = NULL;

CHRCOV *CCOVhead = NULL;
CHRCOV *CCOVcurr = NULL;

RATIOS *rhead = NULL;
RATIOS *rcurr = NULL;

int no_of_samples;
int bamcoverage = 0; //0: quick (index), 1: count all reads

void FreeAllocatedData(void) {
    DestroyBAMstruct(BAMhead);
    DestroyCHROMstruct(CHROMhead, no_of_samples);
    //DestroyRatioStruct(rhead);
}

void ComputeCoverageChIPpeak(CMDINPUT *cmd) {
    PEAK *head = NULL;
    BAMFILES *curr = NULL;
    char *ofile = NULL;
    int ofile_len = 0;

    CHROMhead = ImportChromosomeDataFromBAM(cmd->bamfiles->name, cmd->no_of_samples, cmd->threads);

    if (cmd->libtype == -1) {
        fprintf(stderr, "Detecting library type\n");
        cmd->libtype = DetectLibraryType(cmd->bamfiles);

        if (cmd->libtype == 0) {
            fprintf(stderr, "\tLibrary seems single-end\n");

            if (cmd->fragment_count_mode == 1 && cmd->fragment_size == 0) {
                fprintf(stderr, "ERROR: fragment mode counting is enable, library is single-end, but fragment size is set to 0.");
                fprintf(stderr, "WARNING: Please re-run program without enabling fragment-counting mode, or set fragment size");
                PrintMultiCovMessage(cmd->argv[0]);
                return;
            }
        }
        else
            fprintf(stderr, "\tLibrary seems paired-end\n");
    }

    if (cmd->blacklist_file)
        BlacklistChromosomeFiles(CHROMhead, cmd->blacklist_file);

    PrintBlacklistedChromosomes(CHROMhead, cmd->no_of_samples);

    if (cmd->genome_coverage == 1) {
        MultiGenomeReadCoverage(cmd, CHROMhead);
    }
    else {
        fprintf(stderr, "\nComputing coverage from the idx of BAM files\n");
        GetChromosomeCoveragesIDX(CHROMhead, cmd->bamfiles);
        GetGenomeCoveragesIDX(CHROMhead, cmd->bamfiles);
    }

    if (cmd->blacklist_bed) {
        fprintf(stderr, "Subtracting reads from blaklist BED file ( %s )\n", cmd->blacklist_bed);
        SubtractBlacklistedBEDS(cmd->blacklist_bed, CHROMhead, cmd->bamfiles, cmd->libtype);
    }

    ComputeSamplescales(cmd->bamfiles, CHROMhead, 1);
    head = ReadBED(cmd->bedfile, cmd->threads);
    AllocateReadCovs(head, cmd->no_of_samples);

    curr = cmd->bamfiles;

    while (curr != NULL) {
        fprintf(stderr, "\nSample: %s\n", curr->shortname);
        fprintf(stderr, "\tTotal no. of reads: %d\n", curr->read_coverage);
        fprintf(stderr, "\tLibrary size scale: %.2f\n", curr->scale);
        curr = curr->next;
    }

    fprintf(stderr, "\nProcessing BAM files\n");
    MultiCoverage(cmd->bamfiles, head, cmd);

    if (cmd->outdir)
        ofile_len = strlen(cmd->outdir);

    if (cmd->outprefix)
        ofile_len += strlen(cmd->outprefix);

    ofile_len += strlen("raw_coverages.tsv") + 50;
    ofile = (char *) calloc(ofile_len, sizeof (char));

    if (cmd->outdir) {
        strcat(ofile, cmd->outdir);
        strcat(ofile, "/");
    }

    if (cmd->outprefix) {
        strcat(ofile, cmd->outprefix);
        strcat(ofile, ".");
    }

    strcat(ofile, "raw_coverages.tsv");

    WriteMultiCovsRaw(cmd->bamfiles, head, cmd->no_of_samples, ofile);

    if (ofile)
        free(ofile);

    ofile_len += strlen("TPM_normalized_coverages.tsv" + 1);
    ofile = (char *) calloc(ofile_len + 1, sizeof (char));

    if (cmd->outdir) {
        strcat(ofile, cmd->outdir);
        strcat(ofile, "/");
    }

    if (cmd->outprefix) {
        strcat(ofile, cmd->outprefix);
        strcat(ofile, ".");
    }

    strcat(ofile, "TPM_normalized_coverages.tsv");

    CalculateTPM(cmd->bamfiles, head);
    WriteMultiCovsNormalized(cmd->bamfiles, head, cmd->no_of_samples, ofile);

    if (ofile)
        free(ofile);

    ofile_len += strlen("FPKM_normalized_coverages.tsv" + 1);
    ofile = (char *) calloc(ofile_len + 1, sizeof (char));

    if (cmd->outdir) {
        strcat(ofile, cmd->outdir);
        strcat(ofile, "/");
    }

    if (cmd->outprefix) {
        strcat(ofile, cmd->outprefix);
        strcat(ofile, ".");
    }

    strcat(ofile, "FPKM_normalized_coverages.tsv");

    CalculateFPKM(cmd->bamfiles, head);
    WriteMultiCovsNormalized(cmd->bamfiles, head, cmd->no_of_samples, ofile);

    if (ofile)
        free(ofile);

    ofile_len += strlen("Library_normalized_coverages.tsv" + 1);
    ofile = (char *) calloc(ofile_len + 1, sizeof (char));

    if (cmd->outdir) {
        strcat(ofile, cmd->outdir);
        strcat(ofile, "/");
    }

    if (cmd->outprefix) {
        strcat(ofile, cmd->outprefix);
        strcat(ofile, ".");
    }

    strcat(ofile, "Library_normalized_coverages.tsv");

    CalculateLibScaled(cmd->bamfiles, head);
    WriteMultiCovsNormalized(cmd->bamfiles, head, cmd->no_of_samples, ofile);

    if (ofile)
        free(ofile);

    DeleteBEDs(head);
}

void NormalizeBAMSrna(CMDINPUT *cmd) {
    BAMFILES *curr = NULL;

    cmd->fragment_count_mode = 0;
    fprintf(stderr, "Allocating BINS of size %d for chromosomes\n", cmd->binSize);
    CHROMhead = ImportChromosomeDataFromBAM(cmd->bamfiles->name, cmd->no_of_samples, cmd->threads);

    if (cmd->blacklist_file)
        BlacklistChromosomeFiles(CHROMhead, cmd->blacklist_file);
    
    CHROMhead = ComputeBins(CHROMhead, cmd->binSize);
    CHROMhead = AllocateBins(CHROMhead, cmd->no_of_samples);
    
    cmd->chr = CHROMhead;

    if (cmd->libtype == -1) {
        fprintf(stderr, "Detecting library type\n");
        cmd->libtype = DetectLibraryType(cmd->bamfiles);

        if (cmd->libtype == 0) {
            fprintf(stderr, "\tLibrary seems single-end\n");

            if (cmd->fragment_count_mode == 1 && cmd->fragment_size == 0) {
                fprintf(stderr, "ERROR: fragment mode counting is enable, library is single-end, but fragment size is set to 0.");
                fprintf(stderr, "WARNING: Please re-run program without enabling fragment-counting mode, or set fragment size");
                PrintMultiCovMessage(cmd->argv[0]);
                return;
            }
        } else
            fprintf(stderr, "\tLibrary seems paired-end\n");
    }

    //fprintf(stderr, "\nComputing coverage from the idx of BAM files\n");
    //GetChromosomeCoveragesIDX(CHROMhead, cmd->bamfiles);
    if(strcmp(cmd->scale, INPUTS_CUSTOM) != 0) {
        cmd->bamfiles->scale = 1;
        cmd->bamfiles->genome_scale = 1;
    }
    
    if (cmd->genome_coverage > 0 && strcmp(cmd->scale, INPUTS_NO) != 0 && strcmp(cmd->scale, INPUTS_CUSTOM) != 0) {
        fprintf(stderr, "\nComputing coverage from BAM file\n");
        MultiGenomeCoverage(cmd, CHROMhead);
        ComputeSamplescales(cmd->bamfiles, CHROMhead, 1);
        ScaleGenomeCoverage(cmd->bamfiles, CHROMhead);
    }
    
    curr = cmd->bamfiles;

    while (curr != NULL) {
        fprintf(stderr, "\nSample: %s\n", curr->shortname);
        fprintf(stderr, "\tTotal no. of reads: %d\n", curr->read_coverage);
        fprintf(stderr, "\tLibrary size scale: %.2f\n", curr->scale);
        fprintf(stderr, "\tTotal number of filtered reads: %d\n", curr->filtered_reads);
        fprintf(stderr, "\tBases sequenced: %f\n", curr->base_coverage);
        fprintf(stderr, "\tGenome size: %.f\n", CalculateGenomeSize(CHROMhead));
        fprintf(stderr, "\tGenome scale: %f\n", curr->genome_scale);
        
        if(strcmp(cmd->scale, INPUTS_SMALLEST) == 0)
            curr->genome_scale = curr->scale;
        
        curr = curr->next;
    }

    if(cmd->strandsplit == 0) {
        fprintf(stderr, "\nCreating coverage track for: %s\n", cmd->bamfiles->shortname);
        GetGenomeCoverageRNA(cmd, CHROMhead, returnRNAfilename(cmd));
    }
    else {
        fprintf(stderr, "\nCreating positive coverage track for: %s\n", cmd->bamfiles->shortname);
        cmd->strand = 1;
        GetGenomeCoverageRNA(cmd, CHROMhead, returnRNAfilename(cmd));
        
        fprintf(stderr, "\nCreating negative coverage track for: %s\n", cmd->bamfiles->shortname);
        cmd->strand = -1;
        GetGenomeCoverageRNA(cmd, CHROMhead, returnRNAfilename(cmd));
    }
}

void NormalizeBAMS(CMDINPUT *cmd) {
    BAMFILES *curr = NULL;

    cmd->fragment_count_mode = 0;
    fprintf(stderr, "Allocating BINS of size %d for chromosomes\n", cmd->binSize);
    CHROMhead = ImportChromosomeDataFromBAM(cmd->bamfiles->name, cmd->no_of_samples, cmd->threads);

    if (cmd->blacklist_file)
        BlacklistChromosomeFiles(CHROMhead, cmd->blacklist_file);

    CHROMhead = ComputeBins(CHROMhead, cmd->binSize);
    CHROMhead = AllocateBins(CHROMhead, cmd->no_of_samples);
    cmd->chr = CHROMhead;

    if (cmd->libtype == -1) {
        fprintf(stderr, "Detecting library type\n");
        cmd->libtype = DetectLibraryType(cmd->bamfiles);

        if (cmd->libtype == 0) {
            fprintf(stderr, "\tLibrary seems single-end\n");

            if (cmd->fragment_count_mode == 1 && cmd->fragment_size == 0) {
                fprintf(stderr, "ERROR: fragment mode counting is enable, library is single-end, but fragment size is set to 0.");
                fprintf(stderr, "WARNING: Please re-run program without enabling fragment-counting mode, or set fragment size");
                PrintMultiCovMessage(cmd->argv[0]);
                return;
            }
        } else
            fprintf(stderr, "\tLibrary seems paired-end\n");
    }

    fprintf(stderr, "\nComputing coverage from the idx of BAM files\n");
    GetChromosomeCoveragesIDX(CHROMhead, cmd->bamfiles);
    MultiGenomeBaseCoverage(cmd, CHROMhead);
        
    if(cmd->strandsplit == 1 && strcmp(cmd->scale, INPUTS_CUSTOM) != 0) {
        if(strcmp(cmd->operation, INPUTS_RSTRRNA) == 0 || strcmp(cmd->operation, INPUTS_STRRNA) == 0) {
            cmd->bamfiles->base_coverage = cmd->bamfiles->base_coverage + cmd->bamfiles->next->base_coverage;
            cmd->bamfiles->next->base_coverage = cmd->bamfiles->base_coverage;
        }
    }
    
    if(strcmp(cmd->scale, INPUTS_CUSTOM) != 0)
        ScaleGenomeCoverage(cmd->bamfiles, CHROMhead);

    if (cmd->genome_coverage == 0 && strcmp(cmd->scale, INPUTS_CUSTOM) != 0) {
        GetGenomeCoveragesIDX(CHROMhead, cmd->bamfiles);
        ComputeSamplescales(cmd->bamfiles, CHROMhead, 1);
    }

    curr = cmd->bamfiles;

    while (curr != NULL) {
        fprintf(stderr, "\nSample: %s\n", curr->shortname);
        fprintf(stderr, "\tTotal no. of reads: %d\n", curr->read_coverage);
        fprintf(stderr, "\tLibrary size scale: %.2f\n", curr->scale);
        fprintf(stderr, "\tTotal number of filtered reads: %d\n", curr->filtered_reads);
        fprintf(stderr, "\tBases sequenced: %f\n", curr->base_coverage);
        fprintf(stderr, "\tGenome size: %.f\n", CalculateGenomeSize(CHROMhead));
        fprintf(stderr, "\tGenome scale: %f\n", curr->genome_scale);
        
        if(strcmp(cmd->scale, INPUTS_SMALLEST) == 0)
            curr->genome_scale = curr->scale;
        
        curr = curr->next;
    }

    if (strcmp(cmd->scale, "no") != 0) {
        fprintf(stderr, "\nScaling sample(s)\n");
        MultiGenomeScaler(cmd, CHROMhead);
    }

    if (cmd->smoothBin > 0) {
        if (cmd->tracksmooth == 0 || cmd->tracksmooth == 1) {
            fprintf(stderr, "\nSmoothening signal\n");
            MultiGenomeSmoother(cmd, CHROMhead);
        }
    }
    
    fprintf(stderr, "Printing output BigWig files\n");

    curr = cmd->bamfiles;
    if(cmd->strandsplit == 1)
        cmd->strand = -1;
    
    while (curr != NULL) {
        PrintScaledBigWig(cmd, curr, NULL);
        curr = curr->next;
        
        if(cmd->strandsplit == 1)
            cmd->strand = 1;
    }

    if (strcmp(cmd->operation, "scaled") != 0 && strcmp(cmd->operation, "unscaled") != 0) {
        fprintf(stderr, "Transforming coverage tracks: %s\n", cmd->operation);
        MultiGenomeTransform(cmd, CHROMhead);
        cmd->strand = 0;
        
        if (cmd->smoothBin > 0) {
            if (cmd->tracksmooth == 0 || cmd->tracksmooth == 2) {
                fprintf(stderr, "\nSmoothening transformed signal\n");
                MultiGenomeSmoother(cmd, CHROMhead);
            }
        }
        
        curr = cmd->bamfiles->next;

        while (curr != NULL) {
            PrintScaledBigWig(cmd, curr, cmd->bamfiles->shortname);
            curr = curr->next;
        }
    }
}

void PrintUsage(char *pname) {
    char *ptr = strrchr(pname, '/');
    ptr = ptr ? ptr + 1 : (char *) pname;

    fprintf(stderr, "\nBAMscale: a tool to quantify peaks, and scale sequencing data\n");
    fprintf(stderr, "Version: %s\n", "v0.0.6");

    fprintf(stderr, "\nUsage: %s <command>\n", ptr);
    fprintf(stderr, "\n\tCommands\tDescription\n");
    fprintf(stderr, "\t========\t===========\n");

    fprintf(stderr, "\t   cov\t\tCalculate coverage of BED coordinates in BAM file(s). Outputs are raw read counts, FPKM and TPM normalized values.\n");
    fprintf(stderr, "\t   scale\tConvert BAM files to BigWigs; scale one or multiple files to genome size or to each other.\n");
}

/*
 * 
 */
int main(int argc, char **argv) {
    CMDINPUT *cmd = NULL;
    int found = 0;

    if (argc > 1) {
        if (strcmp(argv[1], INPUTS_COV) == 0) {
            found++;
            cmd = MultiCovParser(argc, argv);

            if (cmd != NULL) {
                no_of_samples = cmd->no_of_samples;
                
                ComputeCoverageChIPpeak(cmd);                    
            }
        }

        if (strcmp(argv[1], INPUTS_SCALE) == 0) {
            found++;
            cmd = ScaleParser(argc, argv);

            if (cmd != NULL) {
                no_of_samples = cmd->no_of_samples;
           
                if(strcmp(cmd->operation, INPUTS_RSTRRNA) == 0 || strcmp(cmd->operation, INPUTS_RNA) == 0 || strcmp(cmd->operation, INPUTS_STRRNA) == 0)
                    NormalizeBAMSrna(cmd);
                
                else
                   NormalizeBAMS(cmd); 
            }
        }

        if (found == 0)
            PrintUsage(argv[0]);
    }
    else {
        PrintUsage(argv[0]);
    }

    if (cmd){
        DestroyCMDinput(cmd);
        FreeAllocatedData();
        free(cmd);
    }

    return (EXIT_SUCCESS);
}

