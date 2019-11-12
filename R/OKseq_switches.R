#args = commandArgs(trailingOnly=TRUE)

if (!require(argparse)) {
  PrintErrorMessage()
  stop("Please install package names \"rtracklayer\"")
}



library(argparse)
parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", default="test.bw", 
                    help="Input OK-seq bigWig file.")

parser$add_argument("-b", "--binsize", type="integer", default=50000, 
                    help="Bin size [default %(default)s]")

parser$add_argument("-l", "--leftmin", type="double", default=-0.15, 
                    help="Mean value of bin on left side [default %(default)s]")

parser$add_argument("-r", "--rightmin", type="double", default=0.15, 
                    help="Mean value of bin on right side [default %(default)s]")

parser$add_argument("-d", "--diff", type="double", default=0.4, 
                    help="Minimum difference between left and right side means [default %(default)s]")

parser$add_argument("-s", "--slidesize", type="double", default=50000, 
                    help="Minimum difference between left and right side means [default %(default)s]")


args <- parser$parse_args()

min_left_value = args$leftmin
min_right_value = args$rightmin
min_diff = args$diff
windowsize = args$binsize
infile = args$input 
sliding_window = args$slidesize


if(length(args) < 1){
  parser$print_help()
  stop("Not enough arguments. Replication timing bigwig file has to be specified (created with BAMscale).")
}

if(!file.exists(infile)) {
  parser$print_help()
  stop(paste0("Input file specified does not exist: \"", infile, "\""))
}

if (!require(rtracklayer)) {
  parser$print_help()
  stop("Please install package names \"rtracklayer\"")
}

library(rtracklayer)

GetMergedElements = function(fbw, typeseg) {
  fbw = fbw[fbw$segments == typeseg,]
  fbw = reduce(fbw)
  fbw$length = width(fbw)
  fbw = fbw[fbw$length > min_segment_length]
  fbw$segments = typeseg
  return(fbw)
}

outfile = gsub(pattern ="\\.bw|\\.bigWig|\\.bigwig", replacement = "", x = infile)
outfile = paste0(outfile, ".OKseq_switches.bed")

include_zeroes_segment_thresholds = 0

bw = import(infile)

dat_bw = as.data.frame(bw)

chrnames = as.character(unique(dat_bw$seqnames))
chrnames = chrnames[grepl(pattern = "_", chrnames) == F]

df <- data.frame(chr=character(),
                 start=numeric(), 
                 end=numeric(), 
                 stringsAsFactors=FALSE) 

binsize = dat_bw[1, "end"] - dat_bw[1, 'start'] + 1
windowsize = windowsize / binsize
rsize = args$binsize


for(j in 1:length(chrnames)) {
  print(chrnames[j])
  dat = dat_bw[dat_bw$seqnames == chrnames[j],]
  endcoord = dat[nrow(dat), "end"]
  
  nbins = (endcoord - rsize) / sliding_window
  lastpos = -1
  reslist = list()
  
  for(i in 1:nbins) {
    startcoord = as.integer(((i-1)*sliding_window)) - 1
    endcoord = as.integer(startcoord + rsize*3) + 2
    
    tdat = dat[dat$start >= (startcoord-1) & dat$end <= (endcoord + 2),]
    tscore = rep(tdat$score, tdat$width)
    
    if(length(tscore) >= rsize*2) {
      lscore = mean(tscore[0:rsize])
      rscore = mean(tscore[(length(tscore)-rsize):length(tscore)])
      
      if(rscore > lscore) {
        if(rscore - lscore >= min_diff & lscore <= min_left_value & rscore >= min_right_value) {
          reslist = append(reslist, list(chrnames[j], startcoord, endcoord))
          lastpos = endcoord
        }
      }
    }
  }
  
  if(length(reslist) > 0) {
    leftpos = data.frame(matrix(unlist(reslist), nrow=length(reslist) / 3, byrow=T),stringsAsFactors=FALSE)
    colnames(leftpos) = c("chr", "start", "end")
    leftpos$start = as.numeric(leftpos$start)
    leftpos$end = as.numeric(leftpos$end)

    if(nrow(df) == 0) {
      df = leftpos
    } else {
      df = merge(df, leftpos, all = T, no.dup = F)
    }
  }
}

leftpos = df
leftpos$tval = "."
leftpos$nval = "."
leftpos$color = "77,0,0"

rightpos = leftpos
rightpos$start = rightpos$start + 2*rsize
rightpos$end = rightpos$start + rsize
rightpos$color = "19,108,0"
leftpos$end = leftpos$start + rsize
merged = rbind(leftpos, rightpos)

interpos = leftpos
interpos$start = interpos$start + rsize
interpos$end = interpos$start + rsize
interpos$color = "240,230,140"
merged = rbind(merged, interpos)

merged = merged[order(merged$chr, merged$start),]
merged = merged[,c("chr", "start", "end", "tval", "tval", "nval", "start", "end", "color")]
write.table(merged, outfile, quote = F, sep = "\t", col.names = F, row.names = F)

