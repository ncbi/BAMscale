/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Inputs.c
 * Author: pongorls
 * 
 * Created on December 11, 2018, 8:00 AM
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <pthread.h>
#include <getopt.h>
#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "Inputs.h"
#include "Definitions.h"
#include "BAMstructs.h"

CMDINPUT *CreateCMDinput(void) {
    CMDINPUT *ptr = (CMDINPUT *) malloc(sizeof (CMDINPUT));
    ptr->bamfiles = NULL;
    ptr->chr = NULL;
    ptr->bedfile = NULL;
    ptr->fragment_count_mode = 0;
    ptr->fragment_size = 0;
    ptr->removeduplicates = 1;
    ptr->libtype = -1;
    ptr->mapq = 0;
    ptr->no_of_samples = 0;
    ptr->nounproper = 1;
    ptr->outdir = NULL;
    ptr->outprefix = NULL;
    ptr->strand = 0;
    ptr->threads = 1;
    ptr->remove_unmapped_pair = 1;
    ptr->max_insert_size = 2000;
    ptr->min_insert_size = 0;
    ptr->fragment_size_filter = 0;
    ptr->binSize = 5;
    ptr->binSizeChange = 0;
    ptr->smoothBin = 0;
    ptr->smoothBinChange = 0;
    ptr->tracksmooth = 1;
    ptr->blacklist_file = NULL;
    ptr->blacklist_bed = NULL;
    ptr->argv = NULL;
    ptr->argc = 0;
    ptr->genome_coverage = 1;
    ptr->strandsplit = 0;
    ptr->custom_scale = NULL;
    
    ptr->filtDiffChr = 1;
    ptr->filtInsSize = 0;

    ptr->normtype = strdup(INPUTS_BASE);
    ptr->scale = strdup(INPUTS_GENOME);
    ptr->operation = strdup(INPUTS_SCALED);

    return ptr;
}

int ParseCustomScaling(CMDINPUT *cmd, char *scales) {
    char *ptr;
    BAMFILES *curr = cmd->bamfiles;
    
    ptr = strtok(scales, ",");
    
    while(ptr != NULL) {
        if(curr != NULL) {
            curr->scale = atof(ptr);
            curr->genome_scale = curr->scale;
            curr = curr->next;
        }
        
        else {
           return -1; 
        }
        
        ptr = strtok(NULL, "\t");
    }
    
    if(curr != NULL) {
        return 1;
    }
    
    return 0;
}

void DestroyCMDinput(CMDINPUT *ptr) {
    if (ptr->bamfiles)
        DestroyBAMstruct(ptr->bamfiles);

    if (ptr->bedfile)
        free(ptr->bedfile);

    if (ptr->outdir)
        free(ptr->outdir);

    if (ptr->outprefix)
        free(ptr->outprefix);

    if (ptr->blacklist_file)
        free(ptr->blacklist_file);

    if (ptr->blacklist_bed)
        free(ptr->blacklist_bed);

    if (ptr->normtype)
        free(ptr->normtype);

    if (ptr->scale)
        free(ptr->scale);

    if (ptr->operation)
        free(ptr->operation);

}

void PrintScaleMessage(char *pname) {
    CMDINPUT *cmd = CreateCMDinput();
    char *ptr = strrchr(pname, '/');
    ptr = ptr ? ptr + 1 : (char *) pname;

    fprintf(stderr, "\nScale one or multiple BAM files\n");
    fprintf(stderr, "Version: %s\n", "v1.0");
    fprintf(stderr, "\nUsage: %s scale [OPTIONS] --bam <BAM_1> (--bam <BAM_2> ... --bam <BAM_N>)\n", ptr);
    fprintf(stderr, "\nOutput: Coverage tracks in BigWig format (un-scaled, scaled, genome scaled)\n");
    fprintf(stderr, "\nRequired options:\n");
    fprintf(stderr, "\t--bam|-i <file>\t\tInput BAM file. This has to be specified at least two times.\n");

    fprintf(stderr, "\nLibrary options:\n");
    fprintf(stderr, "\t--libtype|-l <str>\tSequencing type to be used. Can be: single, paired, and auto (default: autodetect)\n");
    fprintf(stderr, "\t--frag|-f <flag>\tCompute coverage using fragments instead of reads (default: no)\n");
    fprintf(stderr, "\t--fragsize|-a <int>\tFragment size to be used to extend single-end library reads\n");

    fprintf(stderr, "\nNormalization, scaling and operation type:\n");
    fprintf(stderr, "\t--normtype|-y <str>\tType of normalization. (default: %s)\n", cmd->normtype);
    fprintf(stderr, "\t\t\t\tIf no normalization is needed, set '--scale no' argument, the program will disregard this option.\n");
    fprintf(stderr, "\t\t\t\tOptions: \n\t\t\t\t  1) reads: No. of mapped reads/fragments\n\t\t\t\t  2) base: Sum of per-base coverage of reads/fragments\n");
    fprintf(stderr, "\n\t--scale|-k <str>\tMethod to scale samples together. (default: %s)\n", cmd->scale);
    fprintf(stderr, "\t\t\t\tOptions are: \n\t\t\t\t  1) no: no scaling, just calculate coverage\n\t\t\t\t  2) smallest: scale reads to smallest library (multiple-samples only)\n");
    fprintf(stderr, "\t\t\t\t  3) genome: scale samples to 1x genome coverage (only possible with 'base' normalization type)\n\n");
    fprintf(stderr, "\t\t\t\t  4) custom: scale to custom scaling factor (--factor or -F <float> has to be supplied)\n\n");
    
    fprintf(stderr, "\n\t--factor|-F <float>\tScaling factor(s) when \"--scale custom\" normalization is selected.\n");
    fprintf(stderr, "\t\t\t\t  If multiple samples are specified, scaling factors should be comma (\",\") delimited.\n");
    fprintf(stderr, "\t\t\t\t  example in case of three input BAM files: 0.643,0.45667,1.3.\n");
    
    fprintf(stderr, "\n\t--operation|-r <str>\tOperation to perform when scaling samples. Default: %s\n", cmd->operation);
    fprintf(stderr, "\t\t\t\tOptions are: \n\t\t\t\t  1) scaled: output scaled tracks.\n\t\t\t\t  2) unscaled: do not scale files in any way.\n\t\t\t\t  2) log2: log2 transform against first BAM file.\n");
    fprintf(stderr, "\t\t\t\t  3) ratio: coverage ratio against first BAM file.\n\t\t\t\t  4) subtract: subtract coverage against first BAM file.\n");
    fprintf(stderr, "\t\t\t\t  5) rfd: OK-seq RFD calculation\n");
    fprintf(stderr, "\t\t\t\t  6) endseq: strand-specific coverages\n");
    fprintf(stderr, "\t\t\t\t  6) endseqr: strand-specific coverages (reverse strand score is negative)\n");
    fprintf(stderr, "\t\t\t\t  7) reptime: replication timing mode for two BAM files (binsize: 100bp, smoothen: 500 bins)\n");
    fprintf(stderr, "\n\t\t\t\tShort description of settings:\n");
    fprintf(stderr, "\t\t\t\tendseq: generates scaled coverage tracks of positive/negative strands,\n");
    fprintf(stderr, "\t\t\t\t\tand the log2 ratios\n");
    fprintf(stderr, "\n\t\t\t\tendseqr: generates scaled coverage tracks of positive/negative strands,\n");
    fprintf(stderr, "\t\t\t\t\tthe negative strand coverage will be negative, and the log2 ratios are calculated\n");
    fprintf(stderr, "\n\t\t\t\treptime: generates scaled coverage tracks and log2 ratios of two BAM files,\n");
    fprintf(stderr, "\t\t\t\t\tsetting the binsize to 100bp and smoothening smoothen to 500 bins\n");

    fprintf(stderr, "\n\t-S <flag>\t\tOutput strand-specific normalized tracks. One BAM file can be specified only\n");
    
    fprintf(stderr, "\n\t--binsize|-z <int>\tSize of bins for output bigWig/bedgraph generation (default: %d)\n", cmd->binSize);

    fprintf(stderr, "\nSequencing coverage computation options:\n");
    fprintf(stderr, "\t--seqcov|-e <int>\tCompute sequencing coverage from BAM file. (default: '1', count reads while parsing BAM)\n");
    fprintf(stderr, "\t\t\t\tOptions are: \n\t\t\t\t  1) 0: use reads in index (only if normalization is set to 'reads')\n\t\t\t\t  2) 1: count reads while parsing BAM(s)\n");
    fprintf(stderr, "\t\t\t\tWARNING: this option is only useful when 'reads' are used for normalization\n");
    fprintf(stderr, "\n\t--blacklist|-c <file>\tInput file with list of chromosomes to blacklist during scaling analysis\n");
    fprintf(stderr, "\n\t--bedsubtract|-u <int>\tBED file with regions to subtract when computing coverage for normalization\n");
    fprintf(stderr, "\t\t\t\tThese coordinates should not overlap so reads are not counted multiple times\n");
    fprintf(stderr, "\n\t--smoothen|-j <int>\tSmoothen signal by calculating mean of N bins flanking both sides of each bin (default: 0)\n");
    fprintf(stderr, "\t\t\t\tIf set to '0', the signal is not smoothened. To turn on specify a value greater than '0'.\n");
    fprintf(stderr, "\t\t\t\tFor replication timing, a good value is to smoothen to 100k bases. If binSize is 100bp, this would be '1000'\n");
    fprintf(stderr, "\n\t--tracksmooth|-b <int>\tWhich tracks should be smoothened when performing smoothening (default: \'1\' meaning only binned track).\n");
    fprintf(stderr, "\t\t\t\tOptions are: \n\t\t\t\t  1) 0: Smoothen scaled and transformed tracks (log2, ratio or subtracted)");
    fprintf(stderr, "\n\t\t\t\t  2) 1: Smoothen only the scaled sequencing track\n");
    fprintf(stderr, "\t\t\t\t  3) 2: Smoothen only the transformed (log2, ratio or subtract) track\n");


    fprintf(stderr, "\nMapping options:\n");
    fprintf(stderr, "\t--mapq|-q <int>\t\tMinimum (at least) mapping quality (default: %d)\n", cmd->mapq);
    fprintf(stderr, "\t--keepdup|-d <flag>\tKeep duplicated reads (default: no)\n");
    fprintf(stderr, "\t--noproper|-p <flag>\tDo not filter un-proper alignments (default: filter)\n");
    fprintf(stderr, "\t--unmappair|-m <flag>\tDo not remove reads with unmapped pairs\n");
    fprintf(stderr, "\t--minfrag|-g <int>\tMinimum fragment size for read pairs (default: %d)\n", cmd->min_insert_size);
    fprintf(stderr, "\t--maxfrag|-x <int>\tMaximum fragment size for read pairs (default: %d)\n", cmd->max_insert_size);
    fprintf(stderr, "\t--fragfilt|-w <flag>\tFilter reads based on fragment size (default: no)\n");
    fprintf(stderr, "\t--diffchr|-W <flag>\tKeep reads where read pair aligns to different chromosome (default: no)\n");

    fprintf(stderr, "\nOutput options:\n");
    fprintf(stderr, "\t--outdir|-o <str>\tOutput directory name (default: \'.\')\n");

    fprintf(stderr, "\nPerformance options:\n");
    fprintf(stderr, "\t--threads|-t <int>\tNo. of threads to use (default: %d)\n", cmd->threads);

    if (cmd)
        DestroyCMDinput(cmd);
}

CMDINPUT *ScaleParser(int argc, char **argv) {
    CMDINPUT *cmd = CreateCMDinput();
    char *libstr = NULL;
    int error = 0;
    int c;

    BAMFILES *curr;

    cmd->argc = argc;
    cmd->argv = argv;
    argv = argv + 1;
    argc -= 1;

    while (1) {
        static struct option long_options[] = {
            {"bam", required_argument, 0, 'i'},
            {"libtype", required_argument, 0, 'l'},
            {"frag", no_argument, 0, 'f'},
            {"fragsize", required_argument, 0, 'a'},
            {"normtype", required_argument, 0, 'y'},
            {"scale", required_argument, 0, 'k'},
            {"factor", required_argument, 0, 'F'},
            {"binsize", required_argument, 0, 'z'},
            {"seqcov", required_argument, 0, 'e'},
            {"blacklist", required_argument, 0, 'c'},
            {"bedsubtract", required_argument, 0, 'u'},
            {"operation", required_argument, 0, 'r'},
            {"smoothen", required_argument, 0, 'j'},
            {"tracksmooth", required_argument, 0, 'b'},
            {"mapq", required_argument, 0, 'q'},
            {"keepdup", no_argument, 0, 'd'},
            {"noproper", no_argument, 0, 'p'},
            {"unmappair", no_argument, 0, 'm'},
            {"minfrag", required_argument, 0, 'g'},
            {"maxfrag", required_argument, 0, 'x'},
            {"fragfilt", no_argument, 0, 'w'},
            {"outdir", required_argument, 0, 'o'},
            {"prefix", required_argument, 0, 'n'},
            {"threads", required_argument, 0, 't'},
            {"strandsplit", no_argument, 0, 'S'},
            {"diffchr", no_argument, 0, 'W'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "i:l:fa:F:y:k:z:e:c:u:r:j:b:q:dpmg:x:wo:n:t:SW", long_options, &option_index);
        /* Detect the end of the options. */
        if (c == -1)
            break;
        
        char *tmpstr = NULL;
        
        if(optarg)
            tmpstr = strdup(optarg);
        
        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf("option %s", long_options[option_index].name);

                if (tmpstr)
                    printf(" with arg %s", tmpstr);

                printf("\n");
                break;

            case 'i':
                if (access(tmpstr, F_OK) != -1) {
                    if(CheckIndexFile(tmpstr) ==0) {
                       fprintf(stderr, "ERROR: BAM file index (.bai) does not exist ( %s )\n", tmpstr);
                       error++;
                    }
                    
                    else {                   
                        cmd->bamfiles = AddBAMstruct(tmpstr, cmd->bamfiles);
                        cmd->no_of_samples++;
                    }
                } else {
                    fprintf(stderr, "ERROR: BAM file does not exist ( %s )\n", tmpstr);
                    error++;
                }
                break;

            case 'l':
                if (strcmp(tmpstr, INPUTS_SINGLE) != 0 && strcmp(tmpstr, INPUTS_PAIRED) != 0 && strcmp(tmpstr, INPUTS_AUTO) != 0) {
                    fprintf(stderr, "ERROR: unknown library type ( %s )\n", tmpstr);
                    error++;
                } else {
                    libstr = strdup(tmpstr);
                    if (strcmp(libstr, INPUTS_SINGLE) == 0)
                        cmd->libtype = 0;

                    if (strcmp(libstr, INPUTS_PAIRED) == 0)
                        cmd->libtype = 1;
                }
                break;

            case 'f':
                cmd->fragment_count_mode = 1;
                break;

            case 'a':
                if (isdigit(tmpstr[0]))
                    cmd->fragment_size = atoi(tmpstr);

                else {
                    fprintf(stderr, "ERROR: invalid input for fragment size ( %s )\n", tmpstr);
                    error++;
                }
                break;

            case 'y':
                if (cmd->normtype)
                    free(cmd->normtype);

                cmd->normtype = NULL;

                if (strcmp(tmpstr, INPUTS_READS) == 0) {
                    cmd->normtype = strdup(INPUTS_READS);
                } else if (strcmp(tmpstr, INPUTS_BASE) == 0) {
                    cmd->normtype = strdup(INPUTS_BASE);
                } else {
                    fprintf(stderr, "ERROR: invalid normalization type specified ( %s )\n", tmpstr);
                    error++;
                }
                break;

            case 'k':
                if (cmd->scale)
                    free(cmd->scale);

                cmd->scale = NULL;

                if (strcmp(tmpstr, "no") == 0) {
                    cmd->scale = strdup(INPUTS_NO);
                } else if (strcmp(tmpstr, INPUTS_SMALLEST) == 0) {
                    cmd->scale = strdup(INPUTS_SMALLEST);
                } else if (strcmp(tmpstr, INPUTS_GENOME) == 0) {
                    cmd->scale = strdup(INPUTS_GENOME);
                } else if (strcmp(tmpstr, INPUTS_CUSTOM) == 0) {
                    cmd->scale = strdup(INPUTS_CUSTOM);
                } else {
                    fprintf(stderr, "ERROR: invalid scaling type specified ( %s )\n", tmpstr);
                    error++;
                }
                break;

            case 'F':
                if(cmd->custom_scale == NULL)
                    cmd->custom_scale = strdup(tmpstr);
                
                else {
                   fprintf(stderr, "ERROR: custom scaling factor already specified ( %s )\n", cmd->custom_scale);
                   error++;
                }
                
                break;
                
            case 'z':
                cmd->binSize = atoi(tmpstr);
                cmd->binSizeChange = 1;

                break;

            case 'e':
                if (isdigit(optarg[0])) {
                    cmd->genome_coverage = atoi(tmpstr);

                    if (cmd->genome_coverage != 0 && cmd->genome_coverage != 1) {
                        fprintf(stderr, "ERROR: Sequencing coverage flag set to unknown value ( %s )\n", tmpstr);
                        error++;
                    }
                } else {
                    fprintf(stderr, "ERROR: Sequencing coverage flag set to unknown value ( %s )\n", tmpstr);
                    error++;
                }

                break;

            case 'c':
                if (access(tmpstr, F_OK) != -1)
                    cmd->blacklist_file = strdup(tmpstr);

                else {
                    fprintf(stderr, "ERROR: blacklist file with chromosomes to be removed does not exist ( %s )\n", tmpstr);
                    error++;
                }

                break;

            case 'u':
                if (access(optarg, F_OK) != -1)
                    cmd->blacklist_bed = strdup(tmpstr);
                else {
                    fprintf(stderr, "ERROR: blacklist BED file coordinates to be removed does not exist ( %s )\n", tmpstr);
                    error++;
                }

                break;

            case 'r':
                if (cmd->operation)
                    free(cmd->operation);

                cmd->operation = NULL;

                if (strcmp(tmpstr, INPUTS_SCALED) == 0) {
                    cmd->operation = strdup(INPUTS_SCALED);
                } else if (strcmp(tmpstr, INPUTS_UNSCALED) == 0) {
                    cmd->operation = strdup(INPUTS_UNSCALED);
                } else if (strcmp(tmpstr, INPUTS_LOG2) == 0) {
                    cmd->operation = strdup(INPUTS_LOG2);
                } else if (strcmp(tmpstr, INPUTS_RATIO) == 0) {
                    cmd->operation = strdup(INPUTS_RATIO);
                } else if (strcmp(tmpstr, INPUTS_SUBSTRACT) == 0) {
                    cmd->operation = strdup(INPUTS_SUBSTRACT);
                } else if(strcmp(tmpstr, INPUTS_RFD) == 0) {                    
                    cmd->operation = strdup(INPUTS_RFD);
                    cmd->strandsplit = 1;
                } else if(strcmp(tmpstr, INPUTS_END) == 0) {                    
                    cmd->operation = strdup(INPUTS_END);
                    cmd->strandsplit = 1;
                } else if(strcmp(tmpstr, INPUTS_ENDR) == 0) {                    
                    cmd->operation = strdup(INPUTS_ENDR);
                    cmd->strandsplit = 1;
                } else if(strcmp(tmpstr, INPUTS_REP) == 0) {                    
                    cmd->operation = strdup(INPUTS_REP);
                } else if(strcmp(tmpstr, INPUTS_RSTRRNA) == 0) {                    
                    cmd->operation = strdup(INPUTS_RSTRRNA);
                    cmd->strandsplit = 1;
                } else if(strcmp(tmpstr, INPUTS_RNA) == 0) {                    
                    cmd->operation = strdup(INPUTS_RNA);
                }else if(strcmp(tmpstr, INPUTS_STRRNA) == 0) {                    
                    cmd->operation = strdup(INPUTS_STRRNA);
                    cmd->strandsplit = 1;
                }else {
                    fprintf(stderr, "ERROR: invalid operation type specified ( %s )\n", tmpstr);
                    error++;
                }
                break;

            case 'j':
                cmd->smoothBin = atoi(tmpstr);
                cmd->smoothBinChange = 1;

                if (cmd->smoothBin < 0) {
                    fprintf(stderr, "ERROR: no. of bins for smoothening can't be negative ( %d )\n", cmd->smoothBin);
                    error++;
                }

                break;

            case 'b':
                cmd->tracksmooth = atoi(tmpstr);

                if (cmd->tracksmooth < 0 || cmd->tracksmooth > 2) {
                    fprintf(stderr, "ERROR: Unknown value when selecting which tracks to smoothen ( %d )\n", cmd->tracksmooth);
                    error++;
                }

                break;

            case 'q':
                if (isdigit(tmpstr[0]))
                    cmd->mapq = atoi(tmpstr);

                else {
                    fprintf(stderr, "ERROR: invalid input for mapq ( %s )\n", tmpstr);
                    error++;
                }
                break;

            case 'd':
                cmd->removeduplicates = 0;
                break;

            case 'p':
                cmd->nounproper = 0;
                break;

            case 'm':
                cmd->remove_unmapped_pair = 0;
                break;

            case 'w':
                cmd->fragment_size_filter = 1;
                cmd->filtInsSize = 1;

                break;

            case 'g':
                if (isdigit(tmpstr[0])) {
                    cmd->min_insert_size = atoi(tmpstr);
                } else {
                    fprintf(stderr, "ERROR: expected numeric value for minimum insert size ( %s )\n", tmpstr);
                }

                break;

            case 'x':
                if (isdigit(tmpstr[0])) {
                    cmd->max_insert_size = atoi(tmpstr);
                } else {
                    fprintf(stderr, "ERROR: expected numeric value for maximum insert size ( %s )\n", tmpstr);
                }

                break;

            case 'o':
                cmd->outdir = strdup(tmpstr);
                break;

            case 'n':
                cmd->outprefix = strdup(tmpstr);
                break;

            case 't':
                if (isdigit(tmpstr[0]))
                    cmd->threads = atoi(tmpstr);

                else {
                    fprintf(stderr, "ERROR: invalid input for threads ( %s )\n", tmpstr);
                    error++;
                }
                break;
                
            case 'S':
                cmd->strandsplit = 1;
                break;
                
            case 'W':
                cmd->filtDiffChr = 0;
                break;

            case '?':
                /* getopt_long already printed an error message. */
                error++;
                break;

            default:
                abort();
        }
        if (tmpstr)
            free(tmpstr);
    }
    
    if (error > 0) {
        PrintScaleMessage(argv[0]);
        return NULL;
    }
    
    if(strcmp(cmd->scale, INPUTS_CUSTOM) == 0) {
        if(cmd->custom_scale == NULL) {
            fprintf(stderr, "ERROR: Custom scaling was specified, but no factors were specified (-F <FLOAT>).\n");
            error++;
        }
        
        int custom_scale_return = ParseCustomScaling(cmd, cmd->custom_scale);
        
        if(custom_scale_return == 1) {
            fprintf(stderr, "ERROR: There are fewer specified scaling factors than BAM files.\n");
            error++;
        }
        
        if(custom_scale_return == -1) {
            fprintf(stderr, "ERROR: More scaling factors were specified than BAM files.\n");
            error++;
        }
    } 
    
    curr = cmd->bamfiles;

    if (curr == NULL) {
        fprintf(stderr, "ERROR: No BAM file(s) specified\n");
        error++;
    }
    
    else {
        if(cmd->strandsplit == 1) {
            if(strcmp(cmd->operation, INPUTS_RSTRRNA) == 0 || strcmp(cmd->operation, INPUTS_RNA) == 0 || strcmp(cmd->operation, INPUTS_STRRNA) == 0) {
                if(curr->next) {
                    fprintf(stderr, "ERROR: Only one BAM file can be specified in RNA-seq mode!\n");
                    error++;
                }
            }
            
            else {
                if(curr->next) {
                   fprintf(stderr, "ERROR: Only one BAM file can be specified for strand-splitting mode\n");
                   error++;
                }

                else {
                    cmd->bamfiles = AddBAMstruct(curr->name, cmd->bamfiles);
                    cmd->bamfiles->next->scale = cmd->bamfiles->scale;
                    cmd->bamfiles->next->genome_scale = cmd->bamfiles->scale;
                    cmd->no_of_samples++;
                }

                cmd->strand = 1;
            }
        }
        
        
    }

    if (strcmp(cmd->normtype, INPUTS_READS) == 0) {
        if (strcmp(cmd->scale, INPUTS_GENOME) == 0) {
            
            fprintf(stderr, "ERROR: Invalid selection of scaling \'%s\' with normalization \'%s\'.\n", cmd->scale, cmd->normtype);
            fprintf(stderr, "\t'reads' normalization can be paired with \'no\' scaling or \'smallest\' scaling\n");
            error++;
        }
    }

    if (strcmp(cmd->scale, INPUTS_NO) == 0 && strcmp(cmd->operation, INPUTS_UNSCALED) != 0) {
        if(strcmp(cmd->operation, INPUTS_RSTRRNA) != 0 && strcmp(cmd->operation, INPUTS_RNA) != 0 && strcmp(cmd->operation, INPUTS_STRRNA) != 0) {
            fprintf(stderr, "ERROR: Scaling was turned off, but operation is not set to \'unscaled\'.\n");
            fprintf(stderr, "\tEither scale samples, or change operation to \'unscaled\'\n");
            error++;
        }
    }

    if (strcmp(cmd->scale, INPUTS_SMALLEST) == 0 && cmd->no_of_samples <= 1) {
        fprintf(stderr, "ERROR: Scaling to smallest sample, but only one sample was specified. Specify additional samples.\n");
        error++;
    }

    if (strcmp(cmd->operation, INPUTS_LOG2) == 0 || strcmp(cmd->operation, INPUTS_RATIO) == 0 || strcmp(cmd->operation, INPUTS_SUBSTRACT) == 0) {
        if (cmd->no_of_samples <= 1) {
            fprintf(stderr, "ERROR: Operation was set to \'%s\', but only one sample was specified. Specify additional samples.\n", cmd->operation);
            error++;
        }
    }

    argv = argv - 1;
    argc += 1;

    if (libstr) {
        if (strcmp(libstr, INPUTS_SINGLE) == 0 && cmd->fragment_size == 0 && cmd->fragment_count_mode == 1) {
            fprintf(stderr, "ERROR: Fragment size has to be specified for fragment coverage computation with single-end reads\n");
            error++;
        }
    }

    if (cmd->outdir != NULL) {
        DIR* dir = opendir(cmd->outdir);

        if (dir)
            closedir(dir);

        else if (ENOENT == errno) {
            if (mkdir(cmd->outdir, 0700) != 0) {
                fprintf(stderr, "ERROR: Could not create output directory: %s\n", cmd->outdir);
                error++;
            }
        } else {
            fprintf(stderr, "ERROR: Could not create output directory: %s\n", cmd->outdir);
            error++;
        }
    }

    if (cmd->threads < 1) {
        fprintf(stderr, "ERROR: No. of threads can't be smaller than 1\n");
        error++;
    }

    if (error > 0) {
        PrintScaleMessage(argv[0]);
        return NULL;
    }
    
    if (strcmp(cmd->operation, INPUTS_RFD) == 0) {
        if(cmd->binSizeChange == 0)
            cmd->binSize = 1000;
    }
    
    if (strcmp(cmd->operation, INPUTS_REP) == 0) {
        if(cmd->binSizeChange == 0)
            cmd->binSize = 100;
        
        if(cmd->smoothBinChange == 0)
            cmd->smoothBin = 500;
    }
    
    fprintf(stderr, "Computing coverage of regions\n");
    fprintf(stderr, "\tBAMS:\n");

    curr = cmd->bamfiles;

    while (curr != NULL) {
        fprintf(stderr, "\t\t%s\n", curr->name);
        curr = curr->next;
    }
    
    fprintf(stderr, "\n\tMin. mapping qual: %d\n", cmd->mapq);
    fprintf(stderr, "\tRemove duplicates: %d\n", cmd->removeduplicates);
    fprintf(stderr, "\tOnly proper pairs: %d\n", cmd->nounproper);
    fprintf(stderr, "\tBin size: %d\n", cmd->binSize);
    fprintf(stderr, "\tLibrary type: %s\n", libstr ? libstr : INPUTS_AUTO);
    fprintf(stderr, "\tFilter reads based on fragment size: %d\n", cmd->fragment_size_filter);

    if (strcmp(libstr ? libstr : INPUTS_AUTO, INPUTS_AUTO) == 0 && cmd->fragment_size_filter == 1) {
        fprintf(stderr, "\t\t(These will be the fragment size limits in case of paired end)\n");
    }

    if ((strcmp(libstr ? libstr : INPUTS_AUTO, INPUTS_PAIRED) == 0 || libstr == NULL) && cmd->fragment_size_filter == 1) {
        fprintf(stderr, "\t\tMin insert: %d\n", cmd->min_insert_size);
        fprintf(stderr, "\t\tMin insert: %d\n", cmd->max_insert_size);
    }

    fprintf(stderr, "\tCompute fragment coverage: %d\n", cmd->fragment_count_mode);

    if (libstr) {
        if (strcmp(libstr, INPUTS_SINGLE) == 0 && cmd->fragment_count_mode == 1)
            fprintf(stderr, "\tFragment size for single-end: %d\n", cmd->fragment_size);
    }

    if(strcmp(cmd->operation, INPUTS_RSTRRNA) == 0)
        fprintf(stderr, "\tRunning RNA-seq mode: stranded coverages (negative values for reverse strand)\n");
    
    if(strcmp(cmd->operation, INPUTS_RNA) == 0)
        fprintf(stderr, "\tRunning RNA-seq mode: standard RNA-seq coverage\n");
    
    if(strcmp(cmd->operation, INPUTS_STRRNA) == 0)
        fprintf(stderr, "\tRunning RNA-seq mode: stranded coverages (both tracks with positive values)\n");
    
    fprintf(stderr, "\tNo. of threads: %d\n", cmd->threads);

    if (cmd->outdir != NULL)
        fprintf(stderr, "\tOutputting files to directory: %s\n", cmd->outdir);

    if (cmd->outprefix != NULL)
        fprintf(stderr, "\tPrefix for output files is: %s\n", cmd->outprefix);

    if (cmd->blacklist_file)
        fprintf(stderr, "\tBlacklisting chromosomes from file: %s\n", cmd->blacklist_file);

    if (cmd->normtype)
        fprintf(stderr, "\tNormalization: %s\n", cmd->normtype);

    if (cmd->scale)
        fprintf(stderr, "\tScaling type: %s\n", cmd->scale);

    if (cmd->smoothBin > 0)
        fprintf(stderr, "\tSmoothening signal with no. of bins: %d\n", cmd->smoothBin);

    fprintf(stderr, "\tOperation type: %s\n", cmd->operation);

    fprintf(stderr, "\n");
    
    if (libstr)
        free(libstr);
    return cmd;
}

void PrintMultiCovMessage(char *pname) {
    CMDINPUT *cmd = CreateCMDinput();
    char *ptr = strrchr(pname, '/');
    ptr = ptr ? ptr + 1 : (char *) pname;

    fprintf(stderr, "\nCalculate coverage of BED coordinates in BAM file(s)\n");
    fprintf(stderr, "Version: %s\n", "v1.0");
    fprintf(stderr, "\nUsage: %s cov [OPTIONS] --bed <BEDFILE> --bam <BAM_1> (--bam <BAM_2> ... --bam <BAM_N>)\n", ptr);
    fprintf(stderr, "\nOutput: Coverage tables (un-normalized, library-size normalized, FPKM and TPM)\n");
    fprintf(stderr, "\nRequired options:\n");
    fprintf(stderr, "\t--bed|-b <file>\t\tInput BED file\n");
    fprintf(stderr, "\t--bam|-i <file>\t\tInput BAM file. This can be specified multiple times in case of multiple BAM files\n");

    fprintf(stderr, "\nLibrary options:\n");
    fprintf(stderr, "\t--libtype|-l <str>\tSequencing type to be used. Can be: single, paired, and auto (default: autodetect)\n");
    fprintf(stderr, "\t--frag|-f <flag>\tCompute coverage using fragments instead of reads (default: no)\n");
    fprintf(stderr, "\t--strand|-s <flag>\tReads need to have same orientation of peaks (default: unstranded)\n");
    fprintf(stderr, "\t--rstrand|-r <flag>\tReads need to have reverse orientation of peaks (default: unstranded)\n");

    fprintf(stderr, "\nSequencing coverage computation options:\n");
    fprintf(stderr, "\t--seqcov|-e <int>\tCompute sequencing coverage from BAM file quickly using the index (option \'0\'),\n");
    fprintf(stderr, "\t\t\t\tor count number of reads by parsing entire BAM file (slower, but more accurate; set to \'1\' [default])\n");
    fprintf(stderr, "\n\t--blacklist|-c <file>\tInput file with list of chromosomes to blacklist when computing coverage for normalization\n");
    fprintf(stderr, "\n\t--bedsubtract|-u <int>\tBED file with regions to subtract when computing coverage for normalization\n");
    fprintf(stderr, "\t\t\t\tThese coordinates should not overlap so reads are not counted multiple times\n");

    fprintf(stderr, "\nMapping options:\n");
    fprintf(stderr, "\t--mapq|-q <int>\t\tMinimum (at least) mapping quality (default: %d)\n", cmd->mapq);
    fprintf(stderr, "\t--keepdup|-d <flag>\tKeep duplicated reads (default: no)\n");
    fprintf(stderr, "\t--noproper|-p <flag>\tDo not filter un-proper alignments (default: filter)\n");
    fprintf(stderr, "\t--unmappair|-m <flag>\tDo not remove reads with unmapped pairs\n");
    fprintf(stderr, "\t--minfrag|-g <int>\tMinimum fragment size for read pairs (default: %d)\n", cmd->min_insert_size);
    fprintf(stderr, "\t--maxfrag|-x <int>\tMaximum fragment size for read pairs (default: %d)\n", cmd->max_insert_size);
    fprintf(stderr, "\t--fragfilt|-w <flag>\tFilter reads based on fragment size (default: no)\n");
    fprintf(stderr, "\t--diffchr|-W <flag>\tKeep reads where read pair aligns to different chromosome (default: no)\n");
    
    fprintf(stderr, "\nOutput options:\n");
    fprintf(stderr, "\t--outdir|-o <str>\tOutput directory name (default: \'.\')\n");
    fprintf(stderr, "\t--prefix|-n <str>\tOutput prefix for file names (default: none)\n");

    fprintf(stderr, "\nPerformance options:\n");
    fprintf(stderr, "\t--threads|-t <int>\tNo. of threads to use (default: %d)\n", cmd->threads);

    if (cmd)
        DestroyCMDinput(cmd);
}

CMDINPUT *MultiCovParser(int argc, char **argv) {
    CMDINPUT *cmd = CreateCMDinput();
    char *libstr = NULL;
    int error = 0;
    int c;

    BAMFILES *curr;

    cmd->argc = argc;
    cmd->argv = argv;

    argv = argv + 1;
    argc -= 1;

    while (1) {
        static struct option long_options[] = {
            {"mapq", required_argument, 0, 'q'},
            {"keepdup", no_argument, 0, 'd'},
            {"noproper", no_argument, 0, 'p'},
            {"frag", no_argument, 0, 'f'},
            {"strand", no_argument, 0, 's'},
            {"rstrand", no_argument, 0, 'r'},
            {"bed", required_argument, 0, 'b'},
            {"bam", required_argument, 0, 'i'},
            {"libtype", required_argument, 0, 'l'},
            {"threads", required_argument, 0, 't'},
            {"fragsize", required_argument, 0, 'a'},
            {"outdir", required_argument, 0, 'o'},
            {"prefix", required_argument, 0, 'n'},
            {"seqcov", required_argument, 0, 'e'},
            {"blacklist", required_argument, 0, 'c'},
            {"bedsubtract", required_argument, 0, 'u'},
            {"unmappair", no_argument, 0, 'm'},
            {"minfrag", required_argument, 0, 'g'},
            {"maxfrag", required_argument, 0, 'x'},
            {"fragfilt", no_argument, 0, 'w'},
            {"diffchr", no_argument, 0, 'W'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        
        c = getopt_long(argc, argv, "rsfpdq:l:i:b:t:a:o:n:c:u:mg:x:we:W", long_options, &option_index);
        /* Detect the end of the options. */
        if (c == -1)
            break;
        
        char *tmpstr = NULL;
        
        if(optarg)
            tmpstr = strdup(optarg);

        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf("option %s", long_options[option_index].name);

                if (tmpstr)
                    printf(" with arg %s", tmpstr);

                printf("\n");
                break;

            case 'w':
                cmd->fragment_size_filter = 1;

                break;

            case 'g':
                if (isdigit(tmpstr[0])) {
                    cmd->min_insert_size = atoi(tmpstr);
                } else {
                    fprintf(stderr, "ERROR: expected numeric value for minimum insert size ( %s )\n", tmpstr);
                }

                break;

            case 'x':
                if (isdigit(optarg[0])) {
                    cmd->max_insert_size = atoi(tmpstr);
                } else {
                    fprintf(stderr, "ERROR: expected numeric value for maximum insert size ( %s )\n", tmpstr);
                }

                break;

            case 'u':
                if (access(tmpstr, F_OK) != -1)
                    cmd->blacklist_bed = strdup(tmpstr);

                else {
                    fprintf(stderr, "ERROR: blacklist BED file coordinates to be removed does not exist ( %s )\n", tmpstr);
                    error++;
                }


                break;

            case 'c':
                if (access(optarg, F_OK) != -1)
                    cmd->blacklist_file = strdup(tmpstr);

                else {
                    fprintf(stderr, "ERROR: blacklist file with chromosomes to be removed does not exist ( %s )\n", tmpstr);
                    error++;
                }


                break;

            case 'e':
                if (isdigit(tmpstr[0])) {
                    cmd->genome_coverage = atoi(tmpstr);

                    if (cmd->genome_coverage != 0 && cmd->genome_coverage != 1) {
                        fprintf(stderr, "ERROR: Sequencing coverage flag set to unknown value ( %s )\n", tmpstr);
                        error++;
                    }
                } else {
                    fprintf(stderr, "ERROR: Sequencing coverage flag set to unknown value ( %s )\n", tmpstr);
                    error++;
                }

                break;

            case 'o':
                cmd->outdir = strdup(tmpstr);
                break;

            case 'n':
                cmd->outprefix = strdup(tmpstr);
                break;

            case 'l':
                if (strcmp(tmpstr, INPUTS_SINGLE) != 0 && strcmp(tmpstr, INPUTS_PAIRED) != 0 && strcmp(tmpstr, INPUTS_AUTO) != 0) {
                    fprintf(stderr, "ERROR: unknown library type ( %s )\n", tmpstr);
                    error++;
                } else {
                    libstr = strdup(tmpstr);

                    if (strcmp(libstr, INPUTS_SINGLE) == 0)
                        cmd->libtype = 0;

                    if (strcmp(libstr, INPUTS_PAIRED) == 0)
                        cmd->libtype = 1;
                }
                break;

            case 'a':
                if (isdigit(tmpstr[0]))
                    cmd->fragment_size = atoi(tmpstr);

                else {
                    fprintf(stderr, "ERROR: invalid input for fragment size ( %s )\n", tmpstr);
                    error++;
                }
                break;

            case 't':
                if (isdigit(tmpstr[0]))
                    cmd->threads = atoi(tmpstr);

                else {
                    fprintf(stderr, "ERROR: invalid input for threads ( %s )\n", tmpstr);
                    error++;
                }
                break;

            case 'i':
                if (access(tmpstr, F_OK) != -1) {
                    if(CheckIndexFile(tmpstr) ==0) {
                       fprintf(stderr, "ERROR: BAM file index (.bai) does not exist ( %s )\n", tmpstr); 
                       error++;
                    }
                    
                    else {                   
                        cmd->bamfiles = AddBAMstruct(tmpstr, cmd->bamfiles);
                        cmd->no_of_samples++;
                    }
                } else {
                    fprintf(stderr, "ERROR: BAM file does not exist ( %s )\n", tmpstr);
                    error++;
                }
                break;

            case 'b':
                if (access(tmpstr, F_OK) != -1)
                    cmd->bedfile = strdup(tmpstr);

                else {
                    fprintf(stderr, "ERROR: BED file does not exist ( %s )\n", tmpstr);
                    error++;
                }

                break;

            case 'd':
                cmd->removeduplicates = 0;
                break;

            case 'r':
                cmd->strand = -1;
                break;

            case 's':
                cmd->strand = 1;
                break;

            case 'f':
                cmd->fragment_count_mode = 1;
                break;

            case 'p':
                cmd->nounproper = 0;
                break;

            case 'm':
                cmd->remove_unmapped_pair = 0;
                break;

            case 'q':
                if (isdigit(tmpstr[0]))
                    cmd->mapq = atoi(tmpstr);

                else {
                    fprintf(stderr, "ERROR: invalid input for mapq ( %s )\n", tmpstr);
                    error++;
                }
                break;
                
            case 'W':
                cmd->filtDiffChr = 0;
                break;

            case '?':
                error++;
                break;

            default:
                abort();
        }
        free(tmpstr);
    }
    
    argv = argv - 1;
    argc += 1;
    curr = cmd->bamfiles;
    
    if (curr == NULL) {
        fprintf(stderr, "ERROR: No BAM file(s) specified\n");
        error++;
    }

    if (cmd->bedfile == NULL) {
        fprintf(stderr, "ERROR: No BED file specified\n");
        error++;
    }

    if (libstr) {
        if (strcmp(libstr, INPUTS_SINGLE) == 0 && cmd->fragment_size == 0 && cmd->fragment_count_mode == 1) {
            fprintf(stderr, "ERROR: Fragment size has to be specified for fragment coverage computation with single-end reads\n");
            error++;
        }
    }

    if (cmd->outdir != NULL) {
        DIR* dir = opendir(cmd->outdir);

        if (dir)
            closedir(dir);

        else if (ENOENT == errno) {
            if (mkdir(cmd->outdir, 0700) != 0) {
                fprintf(stderr, "ERROR: Could not create output directory: %s\n", cmd->outdir);
                error++;
            }
        } else {
            fprintf(stderr, "ERROR: Could not create output directory: %s\n", cmd->outdir);
            error++;
        }

    }

    if (cmd->threads < 1) {
        fprintf(stderr, "ERROR: No. of threads can't be smaller than 1\n");
        error++;
    }

    if (error > 0) {
        PrintMultiCovMessage(argv[0]);
        return (EXIT_SUCCESS);
    }

    fprintf(stderr, "Computing coverage of regions\n");
    fprintf(stderr, "\tBED: %s\n", cmd->bedfile);
    fprintf(stderr, "\tBAMS:\n");

    curr = cmd->bamfiles;

    while (curr != NULL) {
        fprintf(stderr, "\t\t%s\n", curr->name);
        curr = curr->next;
    }

    fprintf(stderr, "\n\tMin. mapping qual: %d\n", cmd->mapq);
    fprintf(stderr, "\tRemove duplicates: %d\n", cmd->removeduplicates);
    fprintf(stderr, "\tOnly proper pairs: %d\n", cmd->nounproper);
    fprintf(stderr, "\tStrandness: %d\n", cmd->strand);
    fprintf(stderr, "\tLibrary type: %s\n", libstr ? libstr : INPUTS_AUTO);
    
    fprintf(stderr, "\tFilter based on fragment size: %d\n", cmd->fragment_size_filter);
    if (strcmp(libstr ? libstr : INPUTS_AUTO, INPUTS_AUTO) == 0 && cmd->fragment_size_filter == 1) {
        fprintf(stderr, "\t\t(These will be the fragment size limits in case of paired end)\n");
    }

    if ((strcmp(libstr ? libstr : INPUTS_AUTO, INPUTS_PAIRED) == 0 || libstr == NULL) && cmd->fragment_size_filter == 1) {
        fprintf(stderr, "\t\tMin insert: %d\n", cmd->min_insert_size);
        fprintf(stderr, "\t\tMin insert: %d\n", cmd->max_insert_size);
    }

    fprintf(stderr, "\tCompute fragment coverage: %d\n", cmd->fragment_count_mode);

    if (libstr) {
        if (strcmp(libstr, INPUTS_SINGLE) == 0 && cmd->fragment_count_mode == 1)
            fprintf(stderr, "\tFragment size for single-end: %d\n", cmd->fragment_size);
    }

    fprintf(stderr, "\tNo. of threads: %d\n", cmd->threads);

    if (cmd->outdir != NULL)
        fprintf(stderr, "\tOutputing files to dir: %s\n", cmd->outdir);

    if (cmd->outprefix != NULL)
        fprintf(stderr, "\tPrefix for output files is: %s\n", cmd->outprefix);

    if (cmd->blacklist_file)
        fprintf(stderr, "\tBlacklisting chromosomes from file: %s\n", cmd->blacklist_file);

    fprintf(stderr, "\n");

    return cmd;
}