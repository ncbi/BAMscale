
setClass("heatmapper", 
         slots=list(
           bedfiles="character",
           bednames="character",
           bwfiles = "character",
           bwnames = "character",
           subset_peaks = "numeric",
           extend_peaks = "numeric",
           peak_nbin = "numeric",
           quant_extend_peak = "numeric",
           
           ## Clustering
           cluster = "character",
           k_means_clusters = "numeric",
           cluster_binned_distro = "numeric",
           sample_cluster_binned_distro = "numeric",
           auto_hclust_method = "character",
           
           samples_used_arranging = "numeric",
           z_score_dat_arranging = "numeric",
           bed = "GRanges",
           bed_annotation = "data.frame",
           ann_stat_table = "data.frame",
           peak_quants = "data.frame",
           peak_quants_ordered = "data.frame",
           raw_binmats = "list",
           binmats = "list",
           heatmaps = "list",
           
           heatmap_lineplot_means = "list",
           lineplot_max = "numeric",
           lineplot_min = "numeric",
           heatmap_max = "numeric",
           heatmap_min = "numeric",
           heatmap_palette = "character",
           heatmap_reverse_palette = "numeric",
           heatmap_plots = "list",
           combined_heatmaps = "HeatmapList",
           
           gene_annotation = "list"))

Quantify_peaks_old = function(bed, bwfile, exntend, use_extend) {
  print("Quantifying peaks:")
  peak_quants = as.data.frame(matrix(ncol = length(bwfile) + 1, nrow = length(bed)))
  colnames(peak_quants) = c(as.character(seq(1,length(bwfile))), "mean")
  rownames(peak_quants) = bed$pname
  
  for(i in 1:length(bwfile)) {
    cat("\tQuantifying peaks for:", bwfile[i],"\n")
    dyn.load("fastBWreader.so")
    
    if(use_extend == 1) {
      if(length(exntend) == 1) {
        peakcoord = as.data.frame(bed + exntend)
      } else {
        peakcoord = as.data.frame(bed + exntend[i])
      }
    } else {
      peakcoord = as.data.frame(bed)
    }
    quantPeaks <- function(fname, chr, start, end, binned, nbins) {
      .C("readBW_multi", as.character(fname), as.character(chr), as.integer(start), as.integer(end), as.integer(nbins), as.numeric(binned), as.integer(length(chr))) #, 
    }
    
    binned = matrix(0, nrow = nrow(peakcoord), ncol = 1)
    a = quantPeaks(bwfile[i], peakcoord$seqnames, peakcoord$start, peakcoord$end, binned, 1)
    peak_quants[,i] = a[[6]]
  }
  
  if(length(bwfile) > 1) {
    peak_quants$mean = rowMeans(peak_quants[,1:length(bwfile)])
  } else {
    peak_quants$mean = as.numeric(peak_quants[,1])
  }
  
  peak_quants$id = bed$bedid
  cat("\t", "Done importing peaks","\n")
  return(peak_quants)
}

ScoreGrangesBWmeanMean = function(coords, bw) {
  overs = as.data.frame(findOverlaps(coords, bw))
  overs$score = bw[overs$subjectHits]$score
  overs.dt = data.table(overs)
  oversum = overs.dt[,list(score = sum(score)), by='queryHits']
  meanvals = rep(0, length(coords))
  meanvals[oversum$queryHits] = oversum$score
  meanvals = meanvals / width(coords)
  return (meanvals)
}


Quantify_peaks = function(bed, bwfile, exntend, use_extend) {
  print("Quantifying peaks:")
  peak_quants = as.data.frame(matrix(ncol = length(bwfile) + 1, nrow = length(bed)))
  colnames(peak_quants) = c(as.character(seq(1,length(bwfile))), "mean")
  rownames(peak_quants) = bed$pname
  
  for(i in 1:length(bwfile)) {
    cat("\tQuantifying peaks for:", bwfile[i],"\n")
    if(use_extend == 1) {
      if(length(exntend) == 1) {
        peakcoord = resize(x = bed,width = exntend*2, fix = "center")
      } else {
        peakcoord = resize(x = bed,width = exntend[i]*2, fix = "center")
      }
    } else {
      peakcoord = bed
    }
    
    bw = rtracklayer::import(bwfile[i], which = peakcoord)
    peak_quants[,i] = ScoreGrangesBWmeanMean(coords = peakcoord, bw = bw)
  }
  
  if(length(bwfile) > 1) {
    peak_quants$mean = rowMeans(peak_quants[,1:length(bwfile)])
  } else {
    peak_quants$mean = as.numeric(peak_quants[,1])
  }
  
  peak_quants$id = bed$bedid
  cat("\t", "Done importing peaks","\n")
  return(peak_quants)
}

ReadBEDfile = function(bedfiles, subset_peaks, individual_subsetting, blacklist_bed, overlap_bed) {
  bed = GRanges()
  bed_remove = NULL
  bed_overlap = NULL
  
  if(!is.null(blacklist_bed)) {
    bed_remove = read.table(blacklist_bed, sep = "\t")[,1:3]
    colnames(bed_remove) = c("chr", "start", "end")
    bed_remove = makeGRangesFromDataFrame(bed_remove)
  }
  
  if(!is.null(overlap_bed)) {
    bed_overlap = read.table(overlap_bed, sep = "\t")[,1:3]
    colnames(bed_overlap) = c("chr", "start", "end")
    bed_overlap = makeGRangesFromDataFrame(bed_overlap)
  }
  
  for(i in 1:length(bedfiles)) {
    cat("\tImporting BED:", bedfiles[i],"\n")
    tbed = read.table(bedfiles[i], sep = "\t")[,1:3]
    colnames(tbed) = c("chr", "start", "end")
    tbed$bedid = i
    tbed = makeGRangesFromDataFrame(tbed, keep.extra.columns = T)
    
    if(!is.null(bed_remove)) {
      overs = findOverlaps(query = bed_remove, subject = tbed)
      tbed = tbed[-unique(subjectHits(overs))]
    }
    
    if(!is.null(bed_overlap)) {
      overs = findOverlaps(query = bed_overlap, subject = tbed)
      tbed = tbed[unique(subjectHits(overs))]
    }
    
    
    if(individual_subsetting == 1) {
      if(subset_peaks > 0 & length(tbed) > subset_peaks) {
        subsets = sort(sample(seq(1,length(tbed)), subset_peaks))
        tbed = tbed[subsets]
      }
    }
    
    if(i == 1) {
      bed = tbed
    } else {
      bed = c(bed, tbed)
    }
  }
  
  if(subset_peaks > 0 & length(bed) > subset_peaks & individual_subsetting == 0) {
    subsets = sort(sample(seq(1,length(bed)), subset_peaks))
    bed = bed[subsets]
  }
  
  names(bed) = seq(1,length(bed))
  bed = bed
  return (bed)
}

Subset_Binned_to_common_peaks = function(raw_binmats) {
  qhits = list()
  
  for(i in 1:length(raw_binmats)) {
    if(i == 1) {
      qhits = rownames(raw_binmats[[i]])
    } else {
      qhits = intersect(qhits, rownames(raw_binmats[[i]]))
    }
  }
  
  for(i in 1:length(raw_binmats)) {
    raw_binmats[[i]] = raw_binmats[[i]][qhits,]
  }
  
  return (raw_binmats)
}


Bin_coordinates_old = function(bed, bwfile, nbin, exntend) {
  print("Importing and binning coordinates:")
  raw_binmats = list() 
  
  for(i in 1:length(bwfile)) {
    cat("\tImporting peaks for:", bwfile[i],"\n")
    dyn.load("fastBWreader.so")
    
    if(length(exntend) == 1) {
      peakcoord = as.data.frame(bed + exntend)
    } else {
      peakcoord = as.data.frame(bed + exntend[i])
    }
    
    quantPeaks <- function(fname, chr, start, end, binned, nbin) {
      .C("readBW_multi", as.character(fname), as.character(chr), as.integer(start), as.integer(end), as.integer(nbin), as.numeric(binned), as.integer(length(chr))) #, 
    }
    
    binned = matrix(0, nrow = nrow(peakcoord), ncol = nbin)
    a = quantPeaks(bwfile[i], peakcoord$seqnames, peakcoord$start, peakcoord$end, binned, nbin)
    raw_binmats[[i]] = matrix(a[[6]], nrow = nrow(peakcoord), ncol = nbin, byrow = T)
    rownames(raw_binmats[[i]]) = names(bed)
  }
  
  raw_binmats = Subset_Binned_to_common_peaks(raw_binmats)
  cat("\t", "Done quantifying peaks","\n")
  return (raw_binmats)
}

Bin_coordinates = function(bed, bwfile, nbin, exntend) {
  print("Importing and binning coordinates:")
  raw_binmats = list() 
  
  for(i in 1:length(bwfile)) {
    cat("\tImporting peaks for:", bwfile[i],"\n")
    if(length(exntend) == 1) {
      peakcoord = resize(x = bed,width = exntend*2, fix = "center")
    } else {
      peakcoord = resize(x = bed,width = exntend[i]*2, fix = "center")
    }
    
    bw = rtracklayer::import(bwfile[i], which = peakcoord)
    bw.df = as.data.frame(bw)
    bw.df$seqnames = as.character(bw.df$seqnames)
    
    bw.expanded = data.frame(chr = rep(as.character(bw.df$seqnames), times = bw.df$width),
                             start = unlist(sapply(seq(1,nrow(bw.df)), function(i) {
                               bw.df$start[i]:bw.df$end[i]
                             })),
                             score = rep(x = bw.df$score, times = bw.df$width))
    bw.expanded$end = bw.expanded$start
    bw.expanded = makeGRangesFromDataFrame(bw.expanded, keep.extra.columns = T)
    
    tiled.bws = tile(x = peakcoord, n = nbin)
    tiled.bws.unl = unlist(tiled.bws, use.names = F)
    names(tiled.bws.unl) = seq(1,length(tiled.bws.unl))
    tiled.bws.unl$id = rep(x = names(tiled.bws), each = nbin)
    tiled.bws.unl$bin = rep(x = 1:nbin, times = length(tiled.bws))
    tiled.bws.unl$score = ScoreGrangesBWmeanMean(coords = tiled.bws.unl, bw = bw.expanded)
    tiled.bws.unl.df = as.data.frame(tiled.bws.unl)
    
    raw_binmats[[i]] = as.matrix(reshape2::acast(data = tiled.bws.unl.df, id ~ bin, value.var = "score"))
    #rownames(raw_binmats[[i]]) = names(bed)
  }
  
  raw_binmats = Subset_Binned_to_common_peaks(raw_binmats)
  cat("\t", "Done quantifying peaks","\n")
  return (raw_binmats)
}



ImportHeatMapperData = function(bedfiles = NULL,
                                bednames = NULL,
                                bwfiles = NULL,
                                bwnames = NULL,
                                subset_peaks = 10000,
                                peak_nbin = 100,
                                extend_peaks = 5000,
                                quantify_extended_peak = 1,
                                individual_subsetting = 0,
                                blacklist_bed = NULL,
                                overlap_bed = NULL) {
  errors = 0
  
  if(is.null(bedfiles)) {
    print("ERROR: no bedfile specified")
    errors = errors + 1
  }
  
  if(is.null(bwfiles)) {
    print("ERROR: no bigwig specified")
    errors = errors + 1
  }
  
  if(errors > 0) {
    return (NULL)
  }
  
  if(length(extend_peaks) > 1 & length(extend_peaks) != length(bwfiles)) {
    print ("ERROR: number of peak extensions is not equal to number of bigwig files. Specify either one value (for all) or equal number as bigwig files")
    return (NULL)
  }
  
  obj <- new("heatmapper")
  obj@peak_nbin = peak_nbin
  obj@extend_peaks = extend_peaks
  obj@bedfiles = bedfiles
  obj@quant_extend_peak = quantify_extended_peak
  if(is.null(bednames)) {
    obj@bednames = basename(bedfiles)
  } else {
    if(length(bedfiles) != length(bednames)) {
      print(paste0("ERROR: number of specified bed names (",length(bednames),") is not equal to number of bed files(",length(bedfiles),")"))
      errors = errors + 1
    } else {
      obj@bednames = bednames
    }
  }
  
  obj@bwfiles = bwfiles
  
  if(is.null(bwnames)) {
    obj@bwnames = basename(bwfiles)
  } else {
    if(length(bwfiles) != length(bwnames)) {
      print(paste0("ERROR: number of specified bed names (",length(bwnames),") is not equal to number of bed files(",length(bwfiles),")"))
      errors = errors + 1
    } else {
      obj@bwnames = bwnames
    }
  }
  
  if(errors > 0) {
    return (NULL)
  }
  
  for(i in 1:length(obj@bwfiles)) {
    if(!file.exists(obj@bwfiles[i])) {
      print(paste0("ERROR: file \"", obj@bwfiles[i], "\" does not exist"))
      errors = errors + 1
    } 
  }
  
  for(i in 1:length(obj@bedfiles)) {
    if(!file.exists(obj@bedfiles[i])) {
      print(paste0("ERROR: file \"", obj@bedfiles[i], "\" does not exist"))
      errors = errors + 1
    } 
  }
  
  if(errors > 0) {
    return (NULL)
  }
  
  obj@bed = ReadBEDfile(bedfiles = obj@bedfiles, subset_peaks = subset_peaks, individual_subsetting = individual_subsetting, blacklist_bed = blacklist_bed, overlap_bed = overlap_bed)
  
  if(length(obj@bed) > 0) {
    obj@peak_quants = Quantify_peaks(bed = obj@bed, bwfile = obj@bwfiles, exntend = obj@extend_peaks, use_extend = obj@quant_extend_peak)
    obj@raw_binmats = Bin_coordinates(bed = obj@bed, bwfile = obj@bwfiles, nbin = obj@peak_nbin, exntend = obj@extend_peaks)
  }
  
  return (obj)
}

ReorderCLusteredQuantPeaks = function(peak_quants) {
  clusterids = unique(as.numeric(peak_quants$split))
  clustermeans = rep(0, length(clusterids))
  names(clustermeans) = clusterids
  
  for(i in 1:length(clustermeans)) {
    clustermeans[i]=mean(peak_quants[peak_quants$split == clusterids[i], "mean"])
  }
  
  clustermeans = sort(clustermeans, decreasing = T)
  tmp_peak_quants_reordered = data.frame()
  for(i in 1:length(clustermeans)) {
    tmp = peak_quants[peak_quants$split == names(clustermeans[i]),]
    tmp = tmp[order(tmp$mean, decreasing = T),]
    tmp$split = i
    tmp_peak_quants_reordered = rbind(tmp_peak_quants_reordered, tmp)
  }
  
  return (tmp_peak_quants_reordered)
}

Calculate_line_plot_mean_values = function(bwfile, binmats, peak_quants) {
  cmeans = list()
  heatmap_min = c(rep(0, length(bwfile)))
  for(i in 1:length(bwfile)) {
    cmeans[[i]] = as.data.frame(matrix(ncol = length(unique(peak_quants$split)), nrow = ncol(binmats[[1]])))
    colnames(cmeans[[i]]) = unique(peak_quants$split)
    for(j in 1:length(unique(peak_quants$split))) {
      cmeans[[i]][,j] = colMeans(as.matrix(binmats[[i]][rownames(peak_quants[peak_quants$split == unique(peak_quants$split)[j],]),]), na.rm = T)
    }
  }
  
  return (cmeans)
}

Calc_heatmap_max_values = function(binmats, plot_scale_samples_separate) {
  heatmap_max = c(rep(0, length(binmats)))
  
  for(i in 1:length(binmats)) {
    heatmap_max[i] = max(binmats[[i]])#quantile(binmats[[i]],.95, na.rm=T)[[1]]
  }
  
  if(plot_scale_samples_separate == 0) {
    heatmap_max = c(rep(max(heatmap_max), length(binmats)))
  }
  
  return (heatmap_max)
}

Calc_heatmap_min_values  = function(binmats, plot_scale_samples_separate) {
  heatmap_min = c(rep(0, length(binmats)))
  
  for(i in 1:length(binmats)) {
    heatmap_min[i] = quantile(binmats[[i]],.05, na.rm=T)[[1]]
  }
  
  if(plot_scale_samples_separate == 0) {
    heatmap_min = c(rep(max(heatmap_min), length(binmats)))
  }
  
  return (heatmap_min)
}


PrepBEDheatmapdata = function(obj = NULL,
                              cluster = "no",
                              kmeans_clusters = 3,
                              auto_hclust_method = "ward.D",
                              zscore_data_for_arranging = 0,
                              samples_used_arranging = 0) {
  
  obj@cluster = cluster
  obj@k_means_clusters = kmeans_clusters
  obj@auto_hclust_method = auto_hclust_method
  obj@binmats = obj@raw_binmats
  obj@z_score_dat_arranging = zscore_data_for_arranging
  obj@samples_used_arranging = samples_used_arranging
  obj@heatmap_palette = "RdYlBu"
  obj@heatmap_reverse_palette = 0
  
  if(min(samples_used_arranging) < 0 | max(samples_used_arranging) > length(obj@bwfiles)) {
    print(paste0("ERROR: samples selected in samples_used_arranging: ", samples_used_arranging, " variable are smaller or larger than number of bwfiles (", length(obj@bwfiles), ")"))
  }
  
  if(samples_used_arranging[1] == 0) {
    obj@peak_quants_ordered = obj@peak_quants
    
    if(obj@cluster == "no") {
      obj@bed$split = obj@bed$bedid
      obj@peak_quants_ordered$split = obj@peak_quants$id
    }
    
    if(obj@z_score_dat_arranging == 1) {
      for(i in 1:length(obj@bwfiles)) {
        obj@peak_quants_ordered[,i] = as.numeric(scale(obj@peak_quants_ordered[,i]))
      }
    }
    
    if(length(obj@bwfiles) == 1) {
      obj@peak_quants_ordered$mean = obj@peak_quants_ordered[,1]
    } else {
      obj@peak_quants_ordered$mean = rowMeans(obj@peak_quants_ordered[,1:length(obj@bwfiles)])
    }
  } else {
    obj@peak_quants_ordered = obj@peak_quants[,c(obj@samples_used_arranging, "mean", "id")]
    
    if(obj@cluster == "no") {
      obj@bed$split = obj@bed$bedid
      obj@peak_quants_ordered$split = obj@peak_quants$id
    }

    if(obj@z_score_dat_arranging == 1) {
      for(i in 1:length(obj@samples_used_arranging)) {
        obj@peak_quants_ordered[,i] = as.numeric(scale(obj@peak_quants_ordered[,i]))
      }
    }
    
    if(length(obj@samples_used_arranging) == 1) {
      obj@peak_quants_ordered$mean = obj@peak_quants_ordered[,1]
    } else {
      obj@peak_quants_ordered$mean = rowMeans(obj@peak_quants_ordered[,1:length(obj@samples_used_arranging)])
    }
  }
  
  if(obj@cluster == "no") {
    obj@bed$split = obj@bed$bedid
    obj@peak_quants_ordered$split = obj@peak_quants$id
  } else if(obj@cluster == "yes") {
    if(obj@samples_used_arranging[1] == 0) {
      kmeans_clustering = kmeans(obj@peak_quants_ordered[,1:length(obj@bwfiles)], centers = obj@k_means_clusters, algorithm = "Hartigan-Wong", iter.max = 100)
    } else {
      kmeans_clustering = kmeans(obj@peak_quants_ordered[,1:length(obj@samples_used_arranging)], centers = obj@k_means_clusters, algorithm = "Hartigan-Wong", iter.max = 100)
    }
    
    obj@peak_quants_ordered$split = kmeans_clustering$cluster
    obj@peak_quants_ordered = ReorderCLusteredQuantPeaks(obj@peak_quants_ordered)
  } else if(obj@cluster == "auto") {
    if(obj@samples_used_arranging[1] == 0) {
      peak_dist = dist(obj@peak_quants_ordered[,1:length(obj@bwfiles)]) # calculated distance of genes
      peak_hclus <- hclust(peak_dist, method = obj@auto_hclust_method) # cluster genes
      #split genes into clusters
      obj@peak_quants_ordered$split = cutreeDynamic(peak_hclus, 
                                                    method = "hybrid", 
                                                    distM = as.matrix(peak_dist))
    } else {
      peak_dist = dist(obj@peak_quants_ordered[,1:length(obj@samples_used_arranging)]) # calculated distance of genes
      peak_hclus <- hclust(peak_dist, method = obj@auto_hclust_method) # cluster genes
      #split genes into clusters
      obj@peak_quants_ordered$split = cutreeDynamic(peak_hclus, 
                                                    method = "hybrid", 
                                                    distM = as.matrix(peak_dist))
    }
    
    obj@peak_quants_ordered = ReorderCLusteredQuantPeaks(obj@peak_quants_ordered)
  } else {
    print("ERROR: set cluster variable to: no, yes or auto.")
    return(obj)
  }
  
  obj@peak_quants_ordered = obj@peak_quants_ordered[order(-obj@peak_quants_ordered$split, obj@peak_quants_ordered$mean, decreasing = T),]
  obj@heatmap_lineplot_means = Calculate_line_plot_mean_values(obj@bwfiles, obj@binmats, obj@peak_quants_ordered)
  obj@lineplot_max = rep(0,length(obj@heatmap_lineplot_means))
  obj@lineplot_min = rep(0,length(obj@heatmap_lineplot_means))
  
  for(i in 1:length(obj@heatmap_lineplot_means)) {
    obj@lineplot_max[i] = max(obj@heatmap_lineplot_means[[i]])
    obj@lineplot_min[i] = min(obj@heatmap_lineplot_means[[i]])
  }
  
  obj@heatmap_max = Calc_heatmap_max_values(obj@heatmap_lineplot_means, 1)
  obj@heatmap_min = Calc_heatmap_min_values(obj@heatmap_lineplot_means,1)
  
  #if(same_scale == 0) {
  #  obj@heatmap_max = Calc_heatmap_max_values(obj@heatmap_lineplot_means, 1)
  #  obj@heatmap_min = Calc_heatmap_min_values(obj@heatmap_lineplot_means,1)
  #} else {
  #  obj@lineplot_max = rep(max(obj@lineplot_max),length(obj@heatmap_lineplot_means))
  #  obj@lineplot_min = rep(min(obj@lineplot_min),length(obj@heatmap_lineplot_means))
  #  obj@heatmap_max = Calc_heatmap_max_values(obj@heatmap_lineplot_means, 0)
  #  obj@heatmap_min = Calc_heatmap_min_values(obj@heatmap_lineplot_means,0)
  #}
  
  if(obj@cluster == "no") {
    obj@peak_quants_ordered$split = obj@bednames[obj@peak_quants_ordered$split]
    obj@peak_quants_ordered$split = factor(obj@peak_quants_ordered$split, levels = c(obj@bednames))
  }
  
  return(obj)
}

PrepHeatmapPlots = function(obj = NULL,
                            color_palette = "RdYlBu",
                            reverse_palette = 0,
                            heatmap_border = T,
                            lineplot_max_value = NULL,
                            lineplot_min_value = NULL,
                            heatmap_max_value = NULL,
                            heatmap_min_value = NULL,
                            same_scale = 0,
                            split_colors = NULL,
                            raster_quality = 4,
                            extend_line_plot = 0.05) {
  
  obj@combined_heatmaps = HeatmapList()
  
  if(length(color_palette) > 1) {
    if(length(color_palette) != length(obj@binmats)) {
      print("ERROR: specify one color palette, or same number as bigwig files")
      return (obj)
    }
  }
  
  obj@heatmap_palette = color_palette
  brewer_pals = brewer.pal.info
  
  if(length(color_palette[!color_palette %in% rownames(brewer_pals)]) > 0) {
    print(paste0("ERROR: color palette selected: ", color_palette[!color_palette %in% rownames(brewer_pals)], " is not in list of palettes:"))
    print(brewer_pals)
    return (obj)
  }
  
  if(is.null(heatmap_max_value)) {
    heatmap_max_value = Calc_heatmap_max_values(obj@heatmap_lineplot_means, 1)
  } else {
    if(length(heatmap_max_value) == 1) {
      heatmap_max_value = rep(heatmap_max_value, length(obj@bwfiles))
    } else if(length(heatmap_max_value) == length(obj@bwfiles)) {
      
    } else {
      print(paste0("ERROR: incorrect number of heatmap maximums specified. Either specify one value, or same numbers as bigwig files "))
      return (obj)
    }
  }
  
  if(is.null(heatmap_min_value)) {
    heatmap_min_value = Calc_heatmap_min_values(obj@heatmap_lineplot_means, 1)
  } else {
    if(length(heatmap_min_value) == 1) {
      heatmap_min_value = rep(heatmap_min_value, length(obj@bwfiles))
    } else if(length(heatmap_min_value) == length(obj@bwfiles)) {
      
    } else {
      print(paste0("ERROR: incorrect number of heatmap minimums specified. Either specify one value, or same numbers as bigwig files "))
      return (obj)
    }
  }
  
  if(is.null(lineplot_max_value)) {
    lineplot_max_value = obj@lineplot_max
  } else {
    if(length(lineplot_max_value) == 1) {
      lineplot_max_value = rep(lineplot_max_value, length(obj@bwfiles))
    } else if(length(lineplot_max_value) == length(obj@bwfiles)) {
      
    } else {
      print(paste0("ERROR: incorrect number of lineplot maximums specified. Either specify one value, or same numbers as bigwig files "))
      return (obj)
    }
  }
  
  if(is.null(lineplot_min_value)) {
    lineplot_min_value = obj@lineplot_min
  } else {
    if(length(lineplot_min_value) == 1) {
      lineplot_min_value = rep(lineplot_min_value, length(obj@bwfiles))
    } else if(length(lineplot_min_value) == length(obj@bwfiles)) {
      
    } else {
      print(paste0("ERROR: incorrect number of lineplot minimums specified. Either specify one value, or same numbers as bigwig files "))
      return (obj)
    }
  }
  
  if(same_scale == 1) {
    lineplot_min_value = rep(min(lineplot_min_value), length(obj@bwfiles))
    lineplot_max_value = rep(max(lineplot_max_value), length(obj@bwfiles))
    heatmap_min_value = rep(min(heatmap_min_value), length(obj@bwfiles))
    heatmap_max_value = rep(max(heatmap_max_value), length(obj@bwfiles))
  }
  
  if(is.null(split_colors)) {
    split_colors = 1:length(unique(obj@peak_quants_ordered$split))
  } else {
    if(length(split_colors) != length(unique(obj@peak_quants_ordered$split))) {
      print(paste0("ERROR: Incorrect number of colors specified: ", length(split_colors), ". Heatmap has ", length(unique(obj@peak_quants_ordered$split)), " groups."))
      return (obj)
    }
  }
  
  dev_dimensions = round(dev.size("cm"))
  print(dev_dimensions)
  dev_dimensions[2] = dev_dimensions[2] * 0.6
  dev_dimensions[1] = dev_dimensions[1] * 0.6
  
  heatmap_height_ratio = 14
  heatmap_annotation_ratio = 2
  heatmap_width_ratio = 3
  heatmap_gap = 1
  n_plots = length(obj@bwnames)
  
  width_elsize = dev_dimensions[1] / (heatmap_width_ratio*length(obj@bwnames) + heatmap_gap*(length(obj@bwnames)-1))
  height_elsize = width_elsize * (heatmap_height_ratio + heatmap_annotation_ratio)
  
  if(height_elsize > dev_dimensions[2]) {
    size_fact = dev_dimensions[2] / height_elsize
    width_elsize = width_elsize * size_fact
  }
  
  print(width_elsize)
  
  
  #heatmap_annotation_height = dev_dimensions[2]*0.125
  #heatmap_height = dev_dimensions[2] * 0.8
  #heatmap_width = 
  
  for(i in 1:length(obj@binmats)) {
    mean_values = obj@heatmap_lineplot_means[[i]]
    extend_value = (lineplot_max_value[i] - lineplot_min_value[i])*extend_line_plot
    print(paste0(lineplot_max_value[i], ":", lineplot_max_value[i]+extend_value))
    mean_values[mean_values > lineplot_max_value[i]+extend_value] = lineplot_max_value[i] + extend_value
    mean_values[mean_values < lineplot_min_value[i]-extend_value] = lineplot_min_value[i] - extend_value
    
    ha = HeatmapAnnotation(mean = anno_lines(mean_values, 
                                             ylim = c(lineplot_min_value[i]-extend_value,lineplot_max_value[i]+extend_value),
                                             extend = 0,
                                             height = unit(heatmap_annotation_ratio*width_elsize, "cm"),
                                             #gp = gpar(col = split_colors, fontsize = 20*width_elsize))
                                             gp = gpar(col = split_colors),
                                             axis_param = list(
                                               gp = gpar(fontsize = 10*width_elsize)
                                             )),
                           show_annotation_name = c(mean = F))
    
    plot_palette = obj@heatmap_palette[1]
    
    if(length(obj@heatmap_palette) > 1) {
      plot_palette = obj@heatmap_palette[i]
    }
    
    brewer_pals = brewer.pal.info[1]
    maxcols = brewer_pals[plot_palette, "maxcolors"]
    plot_colors = c(brewer.pal(n = maxcols, name = plot_palette)[1], brewer.pal(n = 5, name = plot_palette))
    
    if(obj@heatmap_reverse_palette == 1) {
      plot_colors = rev(plot_colors)
    }
    
    col_fun = colorRamp2(c(heatmap_min_value[i], 
                           heatmap_min_value[i] + heatmap_max_value[i]*0.1,
                           heatmap_min_value[i] + heatmap_max_value[i]*0.2,
                           heatmap_min_value[i] + heatmap_max_value[i]*0.8,
                           heatmap_min_value[i] + heatmap_max_value[i]*1,
                           heatmap_min_value[i] + heatmap_max_value[i]*1.2), plot_colors)
    
    if(heatmap_border == T) {
      ht_opt(heatmap_border = TRUE)
    }
    
    column_names_plot = rep("", ncol(obj@binmats[[i]]))
    column_names_plot[1] = paste0("-",obj@extend_peaks[1], "bp")
    column_names_plot[ncol(obj@binmats[[i]])] = paste0("+", obj@extend_peaks[1], "bp")
    column_names_plot[round(ncol(obj@binmats[[i]])/2)] = "center"
    
    obj@heatmap_plots[[i]] = Heatmap(obj@binmats[[i]][rownames(obj@peak_quants_ordered),],
                                     name = obj@bwnames[i],
                                     column_title = obj@bwnames[i],
                                     column_title_gp = gpar(fontsize = 15*width_elsize),
                                     show_row_names = F,
                                     show_column_names = T,
                                     column_labels = column_names_plot,
                                     column_names_gp = gpar(fontsize = 10*width_elsize),
                                     cluster_rows = F,
                                     cluster_columns = F,
                                     use_raster = T,
                                     raster_quality = raster_quality,
                                     width = unit(heatmap_width_ratio*width_elsize, "cm"), 
                                     height = unit(heatmap_height_ratio*width_elsize, "cm"),
                                     top_annotation = ha,
                                     col = col_fun,
                                     row_split = obj@peak_quants_ordered$split,
                                     row_title_gp = gpar(col = split_colors, font = 2, fontsize = 15*width_elsize),
                                     heatmap_legend_param = list(title_gp = gpar(fontsize = 8*width_elsize), 
                                                                 labels_gp = gpar(fontsize = 8*width_elsize),
                                                                 border = "black",
                                                                 direction = "horizontal",
                                                                 legend_height = unit(8, "cm"),
                                                                 grid_width = unit(1, "cm")))
    
    obj@combined_heatmaps = obj@combined_heatmaps + obj@heatmap_plots[[i]]
  }
  
  ht_opt(heatmap_border = FALSE)
  
  return (obj)
}

AnnotatePeaksBED = function(obj = NULL, txDb = NULL, 
                            tssRegion = c(-3000,3000),
                            level = "transcript",
                            assignGenomicAnnotation = TRUE,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                          "Downstream", "Intergenic"),
                            annoDb = NULL,
                            addFlankGeneInfo = FALSE,
                            flankDistance = 5000,
                            sameStrand = FALSE,
                            ignoreOverlap = FALSE,
                            ignoreUpstream = FALSE,
                            ignoreDownstream = FALSE,
                            overlap = "TSS") {
  if(is.null(obj)) {
    print("ERROR: no object specified")
    return (NULL)
  }
  
  if(is.null(txDb)) {
    print("ERROR: no txdb object specified for ChIPseeker")
    return (obj)
  }
  
  peakAnno <- annotatePeak(peak = obj@bed, tssRegion=tssRegion,
                           TxDb=txDb,
                           level = level, 
                           assignGenomicAnnotation = assignGenomicAnnotation,
                           genomicAnnotationPriority = genomicAnnotationPriority, 
                           annoDb = annoDb, 
                           addFlankGeneInfo = addFlankGeneInfo, 
                           flankDistance = flankDistance, 
                           sameStrand = sameStrand, 
                           ignoreOverlap = ignoreOverlap, 
                           ignoreUpstream = ignoreUpstream, 
                           ignoreDownstream = ignoreDownstream,
                           overlap = overlap)
  
  obj@bed_annotation = as.data.frame(peakAnno)
  
  obj@bed_annotation$annotation_simple = NA
  obj@bed_annotation$annotation_simple = sapply(seq(1,nrow(obj@bed_annotation)), function(i) {
    strsplit(x = obj@bed_annotation$annotation[i], split = " ")[[1]][1]
  })
  
  tmp_annotation = obj@bed_annotation
  tmp_annotation$coord = sapply(seq(1,nrow(tmp_annotation)), function(i) {
    paste0(tmp_annotation$seqnames[i], ":", tmp_annotation$start[i], "-", tmp_annotation$end[i])
  })
  rownames(tmp_annotation) = tmp_annotation$coord
  
  tmp_bed = as.data.frame(obj@bed)
  tmp_bed$name = rownames(tmp_bed)
  tmp_bed$coord = sapply(seq(1,nrow(tmp_bed)), function(i) {
    paste0(tmp_bed$seqnames[i], ":", tmp_bed$start[i], "-", tmp_bed$end[i])
  })
  rownames(tmp_bed) = tmp_bed$coord
  common_coords = intersect(rownames(tmp_annotation), rownames(tmp_bed))
  
  tmp_annotation = tmp_annotation[common_coords,]
  tmp_bed = tmp_bed[common_coords,]
  rownames(tmp_annotation) = tmp_bed$name
  obj@bed_annotation = tmp_annotation
  
  return(obj)
}

PlotAnnotationBarChart = function(obj = NULL, percentage_plot = FALSE) {
  tmpdat = obj@bed_annotation
  pq = obj@peak_quants_ordered
  
  common_feats = intersect(rownames(tmpdat), rownames(pq))
  tmpdat = tmpdat[common_feats,]
  pq = pq[common_feats,]
  tmpdat$split = pq$split
  pcounts = table(pq$split)
  obj@ann_stat_table = as.data.frame(table(tmpdat$split, tmpdat$annotation_simple))
  colnames(obj@ann_stat_table) = c("peak_split", "annotation", "count")
  obj@ann_stat_table$annotation = factor(obj@ann_stat_table$annotation, levels = c("Promoter", "5'", "Exon", "Intron", "3'", "Downstream", "Distal"))
  obj@ann_stat_table$percent = sapply(seq(1,nrow(obj@ann_stat_table)), function(i) {
    obj@ann_stat_table$count[i] * 100 / pcounts[obj@ann_stat_table$peak_split[i]][[1]]
  })
  
  if(percentage_plot == FALSE) {
    p = ggplot(obj@ann_stat_table, aes(x=annotation, fill = peak_split, y = count)) +
      geom_bar(position="dodge", stat="identity") +
      ylab("No. of peaks") +
      xlab("Annotation feature") +
      theme_bw()
  } else {
    p = ggplot(obj@ann_stat_table, aes(x=annotation, fill = peak_split, y = percent)) +
      geom_bar(position="dodge", stat="identity") +
      ylab("Percent of peaks") +
      xlab("Annotation feature") +
      theme_bw()
  }
  
  plot(p)
  
  return (obj)
}

PlotXY = function(obj = NULL,
                  sample_x = 1,
                  sample_y = 2,
                  same_scale = 0,
                  x_limit = NULL,
                  y_limit = NULL,
                  subset_cluster = NULL,
                  percentile = .99,
                  show_lm=T,
                  show_stat=T,
                  point_size = 0.2,
                  cor_test_type = "spearman") {
  
  if(is.character(sample_x)) {
    sample_x = seq(1,length(obj@bwnames))[obj@bwnames == sample_x]
  }
  
  if(is.character(sample_y)) {
    sample_y = seq(1,length(obj@bwnames))[obj@bwnames == sample_y]
  }
  
  plotdf = obj@peak_quants_ordered[,c(sample_x, sample_y)]
  
  if(is.null(x_limit)) {
    x_limit = c(quantile(plotdf[,1], 1-percentile)[[1]], quantile(plotdf[,1], percentile)[[1]])
  }
  
  if(is.null(y_limit)) {
    y_limit = c(quantile(plotdf[,2], 1-percentile)[[1]], quantile(plotdf[,2], percentile)[[1]])
  }
  
  plotdf$split = obj@peak_quants_ordered$split
  
  if(!is.null(subset_cluster)) {
    plotdf = plotdf[plotdf$split %in% subset_cluster,]
  }
  
  colnames(plotdf) = c("x", "y", "split")
  p = ggplot(plotdf, aes(x = x, y = y, color = split)) +
    geom_point(size = point_size) +
    ylim(y_limit) +
    xlim(x_limit) +
    xlab(obj@bwnames[sample_x]) +
    ylab(obj@bwnames[sample_y]) +
    theme_bw()
  
  if(show_stat == T) {
    p = p +stat_cor(method=cor_test_type)
  }
  
  if(show_lm == T) {
    p = p + geom_smooth(method="lm", se=FALSE)
  }
  
  p
  
}