library(rtracklayer)
library(ggplot2)
library(ComplexHeatmap)
require(reshape2)
library(data.table)
library(circlize)
library(RColorBrewer)
source("heatmap_functions_clean.R")

# If you use this script, please cite BAMscale, rtracklayer and ComplexHeatmap!!!
# Disclaimer/FYI: this is still under development!

#Author: Lorinc Pongor (pongorlorinc@gmail.com)

## INPUTS ###

plot_pdf_name = "Heatmap_plot.pdf"

# downsample peaks to these many coordinates (separately each peak)
subset_peaks = 0 # set to 0 to stop
extend = 5000 #extend peaks by these many bases
nbin = 100 #number of bins in heatmap, the default 100 is pretty good

# Number of clusters in k-means clustering.
# Set to 0 to turn it off
k_means_clusters = 3

# Numeric values to set for clustering, or ordering of peaks based on intensity
# Please give numeric values. Eg:
# samples_used_clustering_arranging = c(5,6) for samples 5 and 6
samples_used_clustering_arranging = NULL # set to NULL to use all samples


#Z-score data for k-means clustering or sorting peaks? 0: no, 1: yes
z_score_data_clustering_arranging = 1 

# Scale samples separately?
# Yes: plot_scale_samples_separate = 1
# No: plot_scale_samples_separate = 0
# It is encouraged to set to 1 (No) when comparing bigwigs from different sources.
# If the data is from one experiment (eg. same antibody, same time), and bigwig is normalized, this can be set to 0 
plot_scale_samples_separate = 1

# one or multiple bigwigs
bwfile = c("data/Sample_MCF7_RecQ1_2_020_C_HCVW3BGX2.dd.bam.scaled.bw",
           "data/Sample_MCF7_ER_2_020_C_HCVW3BGX2.dd.bam.scaled.bw",
           "data/FOXA1_ENCFF255FPM.bigWig",
           "data/GATA3_ENCFF477GZL.bigWig",
           "data/H3K27ac_ENCFF411FCW.bigWig",
           "data/H3K4me1_ENCFF983TTS.bigWig",
           "data/H3K4Me3_ENCFF862CKA.bigWig",
           "data/H3K9Me3_ENCFF688REP.bigWig")

#one or multiple peaks
peaks = c("data/MCF7_RecQ1_2_vs_MCF7_IgG_2_MACS2_brdpks.bed")

# Name of bigwigs to be used for the subtitle of each heatmap (column)
# Same number of names have to be specified as the number of bigwig files
# Set to bwnames = NULL to use file name
bwnames = c("RECQ1", "ERa", "FOXA1", "GATA3", "H3K27ac", "H3K4me1", "H3K4Me3", "H3K9Me3")

# Name of peaks to be used for the rowname of each heatmap
# Same number of names have to be specified as the number of peak files
# Set to pnames = NULL to use file name
pnames = NULL

## END OF INPUTS ###


z_score_data = 0 # please don't use this for now

peakds = list()
peakds$peaks = peaks
peakds$bwfile = bwfile
peakds$bwnames = bwnames
peakds$peaknames = pnames
peakds$subset_peaks = subset_peaks
peakds$extend = extend
peakds$nbin = nbin
peakds$k_means_clusters = k_means_clusters
peakds$samples_used_clustering_arranging = samples_used_clustering_arranging
peakds$plot_scale_samples_separate = plot_scale_samples_separate
peakds$z_score_data_clustering_arranging = z_score_data_clustering_arranging
peakds$bed = GRanges()
peakds$peak_quants = data.frame()
peakds$raw_binmats = list()
peakds$binmats = list()
peakds$order = NULL
peakds$clusters = NULL
peakds$errors = 0
peakds$plotlist = list()
peakds$z_score_data = z_score_data
peakds$heatmap_lineplot_means = list()
peakds$max_ann = list()
peakds$min_ann = list()
peakds$heatmap_max = list()
peakds$heatmap_min = list()
source("heatmap_functions_clean.R")
start.time <- Sys.time()
peakds = ImportHeatmapData(peakds)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

if(peakds$errors > 0) {
  stop("Check inputs!!!")
}

peakds = PrepareDataForPlotting(peakds)
peakds = PlotHeatmaps(peakds)
draw(peakds$combined_plot, ht_gap = unit(1, "cm"))

pdf(plot_pdf_name, 
    width = (length(peakds$bwfile) * 4 + (length(peakds$bwfile)-1)) * 0.393701, 
    height = 18 * 0.393701)
draw(peakds$combined_plot, ht_gap = unit(1, "cm"))
dev.off()
