library(GenomicRanges)
library(dynamicTreeCut)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(data.table)

source("Heatmapper_functions.R")


bwfile = c("../../SCLC_NAPY_ChIP/NA_chip/bigwigs/GSM1700639_H889_ASCL1.scaled.bw",
           "../../SCLC_NAPY_ChIP/NA_chip/bigwigs/GSM1700641_H82_NEUROD1.scaled.bw",
           "/Volumes/LMP/ngs/chip/SCLC_cell_lines/bigwigs_hg19_clean/NCI-H1048_POU2F3_rep1.hg19_clean.bam.scaled.bw",
           "/Volumes/LMP/ngs/chip/SCLC_cell_lines/bigwigs_hg19_clean/NCI-H889_H3K27ac_rep1.hg19_clean.bam.scaled.bw",
           "/Volumes/LMP/ngs/chip/SCLC_cell_lines/bigwigs_hg19_clean/NCI-H82_H3K27ac_rep1.hg19_clean.bam.scaled.bw",
           "/Volumes/LMP/ngs/chip/SCLC_cell_lines/bigwigs_hg19_clean/NCI-H1048_H3K27ac_rep1.hg19_clean.bam.scaled.bw",
           "/Volumes/LMP/ngs/chip/SCLC_cell_lines/bigwigs_hg19_clean/DMS114_H3K27ac_rep1.hg19_clean.bam.scaled.bw")

peaks = c("Intervene_results/sets/100_ASCL1.bed",
          "Intervene_results/sets/010_NEUROD1.bed",
          "Intervene_results/sets/001_POU2F3.bed")

bwnames = c("H889_ASCL1", "H82_NEUROD1", "H1048_POU2F3", "H889", "H82", "H1048", "DMS114")
pnames = c("A", "N", "P")

heatobj = ImportHeatMapperData(bedfiles = peaks,
                               bednames = pnames,
                               bwfiles = bwfile, 
                               bwnames = bwnames, 
                               extend_peaks = c(2500),
                               subset_peaks = 1500, 
                               individual_subsetting = 1)

heatobj = PrepBEDheatmapdata(obj = heatobj, cluster = "no")
heatobj = PrepHeatmapPlots(obj = heatobj,
                           color_palette = c("Reds", "Greens", "Blues", "Purples", "Purples", "Purples", "Purples"),
                           #lineplot_max_value = c(20), lineplot_min_value = c(10),
                           split_colors = c("red", "green", "blue"), 
                           raster_quality = 1, same_scale = 0)
draw(heatobj@combined_heatmaps, ht_gap = unit(1, "cm"))

