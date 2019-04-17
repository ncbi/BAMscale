/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BAMstructs.c
 * Author: pongorls
 * 
 * Created on November 28, 2018, 11:58 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "Definitions.h"
#include "BAMstructs.h"

int CheckIndexFile(char *fname) {
    if(fname == NULL)
        return 0;
    
    char *idx = (char *)calloc(strlen(fname) + 5, sizeof(char));
    strcpy(idx, fname);
    strcat(idx, ".bai");
    
    if(access( idx, F_OK ) == -1)
        return 0;
    
    return 1;
}

void DestroyBAMstruct(BAMFILES *head) {
    BAMFILES *curr = head;

    while (head != NULL) {
        curr = head;
        head = head->next;

        if (curr->name)
            free(curr->name);

        if (curr)
            free(curr);
    }

    head = NULL;
    curr = NULL;
}

BAMFILES *AddBAMstruct(char *BAMname, BAMFILES *head) {
    BAMFILES *ptr = (BAMFILES *) calloc (1, sizeof (BAMFILES));
    BAMFILES *curr = head;
    ptr->name = NULL;
    ptr->shortname = NULL;
    ptr->read_coverage = -1;
    ptr->scale = 1;
    ptr->next = NULL;
    ptr->filtered_reads = 0;
    ptr->base_coverage = 0;
    ptr->genome_scale = 0.0;

    ptr->name = strdup(BAMname);
    char *p = strrchr(BAMname, '/');
        
    if(p) {
        p++;
        ptr->shortname = strdup(p);
    }
    
    else {
        ptr->shortname = strdup(BAMname);
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

void PrintBAMstructs(BAMFILES *head) {
    BAMFILES *curr = head;

    while (curr != NULL) {
        printf("File: %s\n\tID:%d\n\tNo. of reads: %d\n\tScale: %f\n", curr->name, curr->id, curr->read_coverage, curr->scale);
        curr = curr->next;
    }
}
