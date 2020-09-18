library(data.table)
library(rtracklayer)

ScoreGrangesBWmean = function(coords, bw) {
  overs = as.data.frame(findOverlaps(coords, bw))
  overs$score = bw[overs$subjectHits]$score
  overs.dt = data.table(overs)
  oversum = overs.dt[,list(score = mean(score)), by='queryHits']
  meanvals = rep(NA, length(coords))
  meanvals[oversum$queryHits] = oversum$score
  return (meanvals)
}

Quantify_peaks = function(bed, bwfile) {
  print("Quantifying peaks:")
  peak_quants = as.data.frame(matrix(ncol = length(bwfile) + 1, nrow = length(bed)))
  colnames(peak_quants) = c(as.character(seq(1,length(bwfile))), "mean")
  rownames(peak_quants) = bed$pname
  
  for(i in 1:length(bwfile)) {
    cat("\tQuantifying peaks for:", bwfile[i],"\n")
    bw = import(bwfile[i], which = bed)
    peak_quants[,i] = ScoreGrangesBWmean(bed, bw)
  }
  
  if(length(bwfile) > 1) {
    peak_quants$mean = rowMeans(peak_quants[,1:length(bwfile)])
  } else {
    peak_quants$mean = as.numeric(peak_quants[,1])
  }
  
  peak_quants$id = bed$id
  cat("\t", "Done importing peaks","\n")
  return(peak_quants)
}

CheckBWfiles = function(peakds) {
  file_not_exist = 0
  if(is.null(length(peakds$bwfile)) | length(peakds$bwfile) < 1) {
    print(paste0("ERROR: no BW files were specified"))
  }
  
  for(i in 1:length(peakds$bwfile)) {
    if(!file.exists(peakds$bwfile[i])) {
      print(paste0("ERROR: file \"", peakds$bwfile[i], "\" does not exist"))
      file_not_exist = 1
    } 
  }
  
  if(file_not_exist == 1) {
    peakds$errors = 1
    return (peakds)
  }
  
  if(is.null(peakds$bwnames)) {
    peakds$bwnames = basename(peakds$bwfile)
  } else {
    if(length(peakds$bwnames) != length(peakds$bwfile)) {
      print(paste0("ERROR: number of specified BW names (",length(peakds$bwnames),") is not equal to number of BW files(",length(peakds$bwfile),")"))
      peakds$errors = 1
    }
  }
  
  if(file_not_exist == 1) {
    peakds$errors = 1
    return (peakds)
  }
  
  return (peakds)
}



CheckBEDfiles = function(peakds) {
  print("Importing BED coordinates")
  file_not_exist = 0
  if(is.null(length(peakds$peaks)) | length(peakds$peaks) < 1) {
    print(paste0("ERROR: no peak files were specified"))
  }
  
  for(i in 1:length(peakds$peaks)) {
    if(!file.exists(peakds$peaks[i])) {
      print(paste0("ERROR: file \"", peakds$peaks[i], "\" does not exist"))
      file_not_exist = 1
    } 
  }
  
  if(file_not_exist == 1) {
    peakds$errors = 1
    return (peakds)
  }
  
  if(is.null(peakds$peaknames)) {
    peakds$peaknames = basename(peakds$peaks)
  } else {
    if(length(peakds$peaknames) != length(peakds$peaks)) {
      print(paste0("ERROR: number of specified peak names (",length(peakds$peaknames),") is not equal to number of peak files(",length(peakds$peaks),")"))
      peakds$errors = 1
    }
  }
  
  if(file_not_exist == 1) {
    peakds$errors = 1
    return (peakds)
  }
  
  for(i in 1:length(peakds$peaks)) {
    cat("\tImporting BED:", peakds$peaks[i],"\n")
    tbed = read.table(peakds$peaks[i], sep = "\t")[,1:3]
    colnames(tbed) = c("chr", "start", "end")
    tbed$id = i
    tbed = makeGRangesFromDataFrame(tbed, keep.extra.columns = T)
    tbed = tbed + extend
    
    if(peakds$subset_peaks > 0 & length(tbed) > peakds$subset_peaks) {
      tbed = tbed[sample(seq(1,length(tbed)), peakds$subset_peaks)]
    }
    
    if(i == 1) {
      bed = tbed
    } else {
      bed = c(bed, tbed)
    }
  }
  
  bed$pname = seq(1, length(bed))
  names(bed) = bed$pname
  peakds$bed = bed
  cat("\t", "Done reading BED files","\n")
  
  return (peakds)
}

BinCoordinates = function(bw, coord, nbins) {
  coord$length = width(coord)
  overs = as.data.frame(findOverlaps(coord, bw))
  overs$width = coord[overs$queryHits]$length
  
  dfcoord = as.data.frame(coord)
  
  dfcoord = dfcoord[overs$queryHits,]
  overs$strand = as.character(dfcoord$strand)
  overs$tss = ifelse(overs$strand == "+", dfcoord$end, dfcoord$start)
  overs$probe_pos = start(bw[overs$subjectHits])
  overs$dist = ifelse(overs$strand == "+", overs$tss - overs$probe_pos, overs$probe_pos - overs$tss)
  overs$bin = round(overs$dist / (overs$width / (nbins-1)))
  overs$bin = overs$bin + 1
  overs$score = bw[overs$subjectHits]$score
  overs.dt = data.table(overs)
  oversum = overs.dt[,list(score = mean(score)), by=c('queryHits', 'bin')]
  
  return (oversum)
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


Bin_peaks = function(bed, bwfile, nbin) {
  print("Importing and binning peaks:")
  raw_binmats = list() 
  
  for(i in 1:length(bwfile)) {
    cat("\tImporting peaks for:", bwfile[i],"\n")
    bw = import(bwfile[i], which = bed)
    df = BinCoordinates(bw, bed, nbin)
    raw_binmats[[i]] = as.data.frame(acast(df,queryHits~bin, value.var = "score"))
    raw_binmats[[i]] = raw_binmats[[i]][,as.character(seq(1:nbin))]
  }
  
  raw_binmats = Subset_Binned_to_common_peaks(raw_binmats)
  cat("\t", "Done quantifying peaks","\n")
  return (raw_binmats)
}

ImportHeatmapData = function(peakds) {
  peakds = CheckBEDfiles(peakds)
  
  if(peakds$errors > 0) {
    return (peakds)
  }
  
  peakds = CheckBWfiles(peakds)
  
  if(peakds$errors > 0) {
    return (peakds)
  }
  
  peakds$peak_quants = Quantify_peaks(peakds$bed, peakds$bwfile)
  peakds$raw_binmats = Bin_peaks(bed = peakds$bed, bwfile = peakds$bwfile, nbin = peakds$nbin)
  peakds$bed = peakds$bed[rownames(peakds$raw_binmats[[1]])]
  
  return (peakds)
}

Scale_binned_matrices = function(raw_binmats, z_score_data) {
  binmats = list()
  
  for(i in 1:length(raw_binmats)) {
    if(z_score_data == 1) {
      binmats[[i]] = as.data.frame(t(scale(t(raw_binmats[[i]]))))
    } else if (z_score_data == -1) {
      binmats[[i]] = as.data.frame(scale(raw_binmats[[i]]))
    } else {
      binmats[[i]] = raw_binmats[[i]]
    }
  }
  
  return (binmats)
}

Arrange_peaks_for_plotting = function(peak_quants, bwfile, k_means_clusters, z_score_data_clustering_arranging, samples_used_clustering_arranging) {
  mean_table = data.frame(mean = peak_quants$mean, id = peak_quants$id)
  tmp_peak_quants = peak_quants # used for ordering temporarily
  
  if(z_score_data_clustering_arranging == 1) {
    for(i in 1:length(bwfile)) {
      tmp_peak_quants[,i] = as.numeric(scale(tmp_peak_quants[,i]))
    }
  }
  
  if(is.null(samples_used_clustering_arranging)) {
    if(length(bwfile) == 1) {
      tmp_peak_quants$mean = as.numeric(tmp_peak_quants[,1])
    } else {
      tmp_peak_quants$mean = rowMeans(tmp_peak_quants[,1:length(bwfile)], na.rm = T)
    }
  } else {
    if(length(samples_used_clustering_arranging) == 1) {
      tmp_peak_quants$mean = as.numeric(tmp_peak_quants[,samples_used_clustering_arranging])
    } else {
      tmp_peak_quants$mean = rowMeans(tmp_peak_quants[,samples_used_clustering_arranging], na.rm = T)
    }
  }
  
  if(k_means_clusters > 0) {
    clusters = NULL
    
    if(length(bwfile) == 1) {
      clusters = kmeans(tmp_peak_quants[,1], centers = k_means_clusters)
    } else if (!is.null(samples_used_clustering_arranging)){
      clusters = kmeans(tmp_peak_quants[,samples_used_clustering_arranging], centers = k_means_clusters)
    } else {
      clusters = kmeans(tmp_peak_quants[,1:length(bwfile)], centers = k_means_clusters)
    }
    
    tmp_peak_quants$split = clusters$cluster
    clusterids = unique(as.numeric(clusters$cluster))
    clustermeans = rep(0, length(clusterids))
    names(clustermeans) = clusterids
    
    for(i in 1:length(clustermeans)) {
      clustermeans[i]=mean(tmp_peak_quants[tmp_peak_quants$split == clusterids[i], "mean"])
    }
    
    clustermeans = sort(clustermeans, decreasing = T)
    peak_quants_reordered = data.frame()
    for(i in 1:length(clustermeans)) {
      tmp = tmp_peak_quants[tmp_peak_quants$split == names(clustermeans[i]),]
      tmp = tmp[order(tmp$mean, decreasing = T),]
      tmp$split = i
      peak_quants_reordered = rbind(peak_quants_reordered, tmp)
    }
    tmp_peak_quants = peak_quants_reordered
    peak_quants = peak_quants[rownames(tmp_peak_quants),]
    peak_quants$split = tmp_peak_quants$split
    
  } else {
    peak_quants$split = peak_quants$id
    tmp_peak_quants = tmp_peak_quants[order(tmp_peak_quants$split, tmp_peak_quants$mean, decreasing = T),]
    peak_quants = peak_quants[rownames(tmp_peak_quants),]
  }
  
  return (peak_quants)
}

Calc_mean_values_for_plots = function(bwfile, binmats, peak_quants) {
  cmeans = list()
  heatmap_min = c(rep(0, length(bwfile)))
  for(i in 1:length(bwfile)) {
    cmeans[[i]] = as.data.frame(matrix(ncol = length(unique(peak_quants$split)), nrow = nbin))
    colnames(cmeans[[i]]) = unique(peak_quants$split)
    for(j in 1:length(unique(peak_quants$split))) {
      cmeans[[i]][,j] = colMeans(as.matrix(binmats[[i]][rownames(peak_quants[peak_quants$split == unique(peak_quants$split)[j],]),]), na.rm = T)
    }
  }
  
  return (cmeans)
}

Calc_max_line_plot = function(cmeans, plot_scale_samples_separate) {
  max_ann = c(rep(0, length(cmeans)))
  
  for(i in 1:length(cmeans)) {
    max_ann[i] = max(cmeans[[i]], na.rm = T)
  }
  
  if(plot_scale_samples_separate == 0) {
    max_ann = c(rep(max(max_ann), length(cmeans)))
  }
  
  return (max_ann)
}

Calc_min_line_plot = function(cmeans, plot_scale_samples_separate) {
  min_ann = c(rep(0, length(cmeans)))
  
  for(i in 1:length(cmeans)) {
    min_ann[i] = min(0, min(cmeans[[i]], na.rm = T))
  }
  
  if(plot_scale_samples_separate == 0) {
    min_ann = c(rep(min(min_ann), length(cmeans)))
  }
  
  return (min_ann)
}

Calc_heatmap_max_values = function(binmats, plot_scale_samples_separate) {
  heatmap_max = c(rep(0, length(binmats)))
  
  for(i in 1:length(binmats)) {
    heatmap_max[i] = quantile(binmats[[i]],.95, na.rm=T)[[1]]
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


PrepareDataForPlotting = function(peakds) {
  peakds$binmats = Scale_binned_matrices(raw_binmats = peakds$raw_binmats, z_score_data = peakds$z_score_data)
  peakds$peak_quants = Arrange_peaks_for_plotting(peak_quants = peakds$peak_quants, 
                                           bwfile = peakds$bwfile,
                                           k_means_clusters = peakds$k_means_clusters, 
                                           z_score_data_clustering_arranging = peakds$z_score_data_clustering_arranging, 
                                           samples_used_clustering_arranging = peakds$samples_used_clustering_arranging)
  
  
  peakds$heatmap_lineplot_means = Calc_mean_values_for_plots(peakds$bwfile, peakds$binmats, peakds$peak_quants)
  peakds$max_ann = Calc_max_line_plot(peakds$heatmap_lineplot_means, peakds$plot_scale_samples_separate)
  peakds$min_ann = Calc_min_line_plot(peakds$heatmap_lineplot_means, peakds$plot_scale_samples_separate)
  peakds$heatmap_max = Calc_heatmap_max_values(binmats = peakds$binmats, peakds$plot_scale_samples_separate)
  peakds$heatmap_min = Calc_heatmap_min_values(binmats = peakds$binmats, peakds$plot_scale_samples_separate)
  
  return (peakds)
}

PlotHeatmaps = function(peakds) {
  peakds$plotlist = list()
  for(i in 1:length(peakds$bwfile)) {
    ptable = peakds$binmats[[i]]
    ha = HeatmapAnnotation(mean = anno_lines(peakds$heatmap_lineplot_means[[i]], 
                                             ylim = c(peakds$min_ann[i],peakds$max_ann[i]), 
                                             height = unit(2, "cm"),
                                             gp = gpar(col = 1:length(unique(peakds$peak_quants$split)))),
                           show_annotation_name = c(mean = FALSE))
    
    col_fun = colorRamp2(c(peakds$heatmap_min[i], 
                           peakds$heatmap_min[i] + peakds$heatmap_max[i]*0.15, 
                           peakds$heatmap_min[i] + peakds$heatmap_max[i]*0.3, 
                           peakds$heatmap_min[i] + peakds$heatmap_max[i]*0.6, 
                           peakds$heatmap_min[i] + peakds$heatmap_max[i], 
                           peakds$heatmap_min[i] + peakds$heatmap_max[i]*1.15), brewer.pal(n = 11, name = "RdYlBu")[c(1,4,6,8,10,11)])
    
    peakds$plotlist[[i]] = Heatmap(ptable[rownames(peakds$peak_quants),],
                            name = peakds$bwnames[i],
                            column_title = peakds$bwnames[i],
                            show_row_names = F,
                            show_column_names = F,
                            cluster_rows = F,
                            cluster_columns = F,
                            use_raster = T,
                            width = unit(3, "cm"), 
                            height = unit(14, "cm"),
                            top_annotation = ha,
                            col = col_fun,
                            row_split = peakds$peak_quants$split,
                            row_title_gp = gpar(col = 1:length(unique(peakds$peak_quants$split)), font = 2))
    
  
    if(i == 1) {
      oplot =  peakds$plotlist[[i]]
    } else {
      oplot = oplot +  peakds$plotlist[[i]]
    }  
  }
  peakds$combined_plot = oplot
  return (peakds)
}