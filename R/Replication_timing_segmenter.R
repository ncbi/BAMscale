args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
  stop("Not enough arguments. Replication timing bigwig file has to be specified (created with BAMscale).")
}

library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)

infile = args[1]
outfile = gsub(pattern ="\\.bw|\\.bigWig|\\.bigwig", replacement = "", x = infile)
outfile = paste0(outfile, ".replication_timings.bed")

options(warn=-1)
GetMergedElements = function(fbw, typeseg) {
  fbw = fbw[fbw$segments == typeseg,]
  fbw = reduce(fbw)
  fbw$length = width(fbw)
  fbw = fbw[fbw$length > min_segment_length]
  fbw$segments = typeseg
  return(fbw)
}

include_zeroes_segment_thresholds = 0
min_segment_length = 5000

bw = import(infile)

if(include_zeroes_segment_thresholds == 0) {
  medval = quantile(bw$score[bw$score != 0], .5, na.rm = T)[[1]]
  upper = quantile(bw$score[bw$score != 0], .75, na.rm = T)[[1]]
  lower = quantile(bw$score[bw$score != 0], .25, na.rm = T)[[1]]
} else {
  medval = quantile(bw$score, .5, na.rm = T)[[1]]
  upper = quantile(bw$score, .75, na.rm = T)[[1]]
  lower = quantile(bw$score, .25, na.rm = T)[[1]]
}


bw$segments = 0
bw$segments = ifelse(bw$score < lower, -2, bw$segments)
bw$segments = ifelse(bw$score < medval & bw$segments == 0, -1, bw$segments)
bw$segments = ifelse(bw$score > medval & bw$segments == 0, 1, bw$segments)
bw$segments = ifelse(bw$score > upper, 2, bw$segments)
bw = resize(bw, width(bw) + 1, fix="start")

segments = GetMergedElements(bw, 2)
segments = c(segments, GetMergedElements(bw, 1))
segments = c(segments, GetMergedElements(bw, -1))
segments = c(segments, GetMergedElements(bw, -2))
segments = sort(segments)
segments = resize(segments, width(segments) - 1, fix="end")

seg_gaps = setdiff(as(seqinfo(segments), "GRanges"), segments)
seg_gaps$length = width(seg_gaps)
seg_gaps$segments = 0
segments = c(segments, seg_gaps)
segments = sort(segments)

if(segments$segments[1] == 0) {
  segments$segments[1] = segments$segments[2]
}

if(segments$segments[length(segments)] == 0) {
  segments$segments[length(segments)] = segments$segments[length(segments)-1]
}

segments$final_seg = sapply(seq(1, length(segments), 1), function(i) {
  if(i > 1 & i < length(segments)) {
    if(segments$segments[i] == 0) {
      if(segments$length[i-1] > segments$length[i+1]) {
        segments$segments[i-1]
      } else {
        segments$segments[i+1]
      }
    } else {
      segments$segments[i]
    }
  } else {
    segments$segments[i]
  }
})

segments = resize(segments, width(segments) + 1, fix="start")
segments$segments = segments$final_seg

segments_final = GetMergedElements(segments, 2)
segments_final = c(segments_final, GetMergedElements(segments, 1))
segments_final = c(segments_final, GetMergedElements(segments, -1))
segments_final = c(segments_final, GetMergedElements(segments, -2))
segments_final = sort(segments_final)
segments_final = resize(segments_final, width(segments_final) - 1, fix="end")

segments_final = as.data.frame(segments_final)
segments_final$width = NULL
segments_final$strand = NULL
segments_final$length = NULL

segments_final$color = "126,0,21"
segments_final$color = ifelse(segments_final$segments == -1, "153,89,31", segments_final$color)
segments_final$color = ifelse(segments_final$segments == 1, "154,161,14", segments_final$color)
segments_final$color = ifelse(segments_final$segments == 2, "20,155,3", segments_final$color)

segments_final$time = "late"
segments_final$time = ifelse(segments_final$segments == -1, "mid-late", segments_final$time)
segments_final$time = ifelse(segments_final$segments == 1, "mid-early", segments_final$time)
segments_final$time = ifelse(segments_final$segments == 2, "early", segments_final$time)

segments_final$tval = 0
segments_final$nval = "."

segments_final = segments_final[,c("seqnames", "start", "end", "time", "tval", "nval", "start", "end", "color")]
write.table(segments_final, outfile, sep = "\t", quote = F, col.names = F, row.names = F)
