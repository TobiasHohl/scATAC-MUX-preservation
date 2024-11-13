
ATAC_QC <- function(seu, save = NULL, filtered = FALSE, plotdir = '', cut_PRFs_low = 3000, cut_PRFs_high = 20000, cut_RIPr = 15, cut_BR = 0.05, cut_NS_upper = 4, cut_NS_lower = 0, cut_TSSenr = 2, region = "chr1-1-2000000", legend.position = 'right', design = 'default', w = 10, h = 10, palette = "Dark2", manual = F){
  # plot QC metrics and include cutoff values
  # metadata  - extracted metadata of seurat object. Object needs to be processed by ComputeQCMetrics.
  # save      - set to one of the following file endings to save plots in plotdir directory: 'eps', 'ps', 'tex', 'pdf', 'jpeg', 'tiff', 'png', 'bmp', 'svg'
  # filtered  - set to TRUE if metadata of filtered seurat is used
  # plotdir   - directory in which the plot png will be saved
  # cut_XYZ   - filter values, thresholds that will be added to the QC plots as horizontal or vertical lines accordingly
  
  metadata <- seu@meta.data
  
  # check for presence of QC metrics
  if(sum(c(grepl('pct_reads_in_peaks', colnames(metadata)),
           grepl('peak_region_fragments', colnames(metadata)),
           grepl('TSS.enrichment', colnames(metadata)),
           grepl('nucleosome_signal', colnames(metadata)),
           grepl('nucleosome_group', colnames(metadata)),
           grepl('high.tss', colnames(metadata)))) != 6) {
    stop('QC values missing. Please first run ComputeQCMetrics on seurat object!')
  }
  
  
  # Cells recovered per sample
  p_cells <- metadata %>% 
    ggplot(aes(x=sample, fill=sample)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("Cells per sample") +
    theme(legend.position = "none", axis.title.x = element_blank())
  # plot ratio of reads in peaks
  p1 <- metadata %>% 
    ggplot(aes(color=sample, x=pct_reads_in_peaks, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    geom_vline(xintercept = cut_RIPr) +
    theme(plot.title = element_text(hjust=0.5, face="bold"), legend.position = legend.position, axis.title.x = element_blank()) +
    ggtitle('FRIPs (%)')
  # fragments in peaks
  if(is.numeric(try(seu$atac_peak_region_fragments, silent = T))) {
    p2 <- metadata %>% 
      ggplot(aes(color=sample, x=atac_peak_region_fragments, fill= sample)) + 
      geom_density(alpha = 0.2) + 
      theme_classic() +
      scale_x_log10() + 
      geom_vline(xintercept = cut_PRFs_low) +
      geom_vline(xintercept = cut_PRFs_high) +
      theme(plot.title = element_text(hjust=0.5, face="bold"), legend.position = legend.position, axis.title.x = element_blank()) +
      ggtitle('Reads in peaks')
  } else {
    p2 <- metadata %>% 
      ggplot(aes(color=sample, x=peak_region_fragments, fill= sample)) + 
      geom_density(alpha = 0.2) + 
      theme_classic() +
      scale_x_log10() + 
      geom_vline(xintercept = cut_PRFs_low) +
      geom_vline(xintercept = cut_PRFs_high) +
      theme(plot.title = element_text(hjust=0.5, face="bold"), legend.position = legend.position, axis.title.x = element_blank()) +
      ggtitle('Reads in peaks')
  }

  # TSS enrichment score per cell
  p3 <- metadata %>% 
    ggplot(aes(color=sample, x=TSS.enrichment, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = cut_TSSenr) +
    theme(plot.title = element_text(hjust=0.5, face="bold"), legend.position = legend.position, axis.title.x = element_blank()) +
    ggtitle('TSS enrichment score')
  # nucleosome signal score
  p5 <- metadata %>% 
    ggplot(aes(color=sample, x=nucleosome_signal, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = cut_NS_upper) +
    geom_vline(xintercept = cut_NS_lower) +
    theme(plot.title = element_text(hjust=0.5, face="bold"), legend.position = legend.position, axis.title.x = element_blank()) +
    ggtitle('Nucleosome signal score')
  # fragment length periodicity
  histo_ns <- FragmentHistogram(object = seu, group.by = "sample", region = region)+ ggtitle('DNA fragment sizes') + theme(plot.title = element_text(hjust=0.5, face="bold"))
  # accessibility signal over all TSS sites
  TSS_groups <- TSSPlot(seu, group.by = 'sample') + NoLegend() + ggtitle('Overall TSS accessibility') + theme(plot.title = element_text(hjust=0.5, face="bold"))
  
  if(manual == TRUE){
    if(length(palette) < length(unique(metadata$sample))){
      stop('Not enough colors given for manual palettes. Make sure to match sample and color amounts!')
    }
    p_cells <- p_cells + scale_fill_manual(values = palette)
    p1 <- p1 + scale_fill_manual(values = palette) + scale_color_manual(values = palette)
    p2 <- p2 + scale_fill_manual(values = palette) + scale_color_manual(values = palette)
    p3 <- p3 + scale_fill_manual(values = palette) + scale_color_manual(values = palette)
    p5 <- p5 + scale_fill_manual(values = palette) + scale_color_manual(values = palette)
    TSS_groups <- TSS_groups + scale_color_manual(values = palette)
    histo_ns <- histo_ns + scale_fill_manual(values = palette)
  } else if (manual == FALSE){
    p_cells <- p_cells + scale_fill_brewer(palette = palette)
    p1 <- p1 + scale_fill_brewer(palette = palette) + scale_color_brewer(palette = palette)
    p2 <- p2 + scale_fill_brewer(palette = palette) + scale_color_brewer(palette = palette)
    p3 <- p3 + scale_fill_brewer(palette = palette) + scale_color_brewer(palette = palette)
    p5 <- p5 + scale_fill_brewer(palette = palette) + scale_color_brewer(palette = palette)
    TSS_groups <- TSS_groups + scale_color_brewer(palette = palette)
    histo_ns <- histo_ns + scale_fill_brewer(palette = palette)
  }
  
  
  if(design == 'default') {
    if(length(unique(seu$sample)) == 1) {
      # create plot layout for patchwork
      design <- '
        12466
        13577
      '
      w <- 12
      h <- 5
    } else {
      # create plot layout for patchwork
      design <- '
        123
        145
        666
        777
      '
      w <- 10
      h <- 10
    }
  }
  
  # combine plots
  plots_QC <- p_cells + p1 + p2 + p3 + p5 + TSS_groups + histo_ns + plot_layout(design = design)
  cat('Plots generated\n')
  
  #print plots into RStudio plots device
  print(plots_QC)
  
  # save plots if save_png is T
  if(!is.null(save)) {
    file_types <- c('eps', 'ps', 'tex', 'pdf', 'jpeg', 'tiff', 'png', 'bmp', 'svg')
    if(!(save %in% file_types)) {
      stop('File extension for saving plots not recognized! Please refer to ggsave to select proper extension.')
    }
    
    # add filtered tag to output file name if applicable
    if(filtered == TRUE) {
      filtered <- '_filtered'
    } else {
      filtered <- ''
    }
    
    ggsave(filename = paste0(plotdir, 'QC-plots', filtered, '.', save),
           plot = plots_QC,
           width = w,
           height = h)
  }
}

ComputeQCMetrics.scATAC <- function(seu, TSS_threshold = 2, nucleosome_threshold = 4, sample_names = NULL, sample_prefixes = NULL) {
  # Compute metrics used for assessing data quality
  # Metrics are:
  # - Nucleosome score
  # - TSS enrichment score
  # - Reads in peaks ratio
  # - FrIPs
  # - Blacklist region fragment ratio
  
  # check for presence of metadata
  if(sum(c(grepl('passed_filters', colnames(seu@meta.data)),
           grepl('peak_region_fragments', colnames(seu@meta.data)),
           grepl('blacklist_region_fragments', colnames(seu@meta.data)))) != 3) {
    stop('Metadata missing! Please check for correct addition of metadata when creating seurat objects!')
  }
  
  # compute nucleosome signal score per cell
  seu <- NucleosomeSignal(object = seu)
  
  # compute TSS enrichment score per cell
  seu <- TSSEnrichment(object = seu, fast = FALSE)
  
  metadata <- seu@meta.data
  # Add cell IDs to metadata
  metadata$cells <- rownames(metadata)
  
  # add blacklist ratio and fraction of reads in peaks
  metadata$pct_reads_in_peaks <- metadata$peak_region_fragments / metadata$passed_filters * 100
  metadata$blacklist_ratio <- metadata$blacklist_region_fragments / metadata$peak_region_fragments
  # group cells in high/low TSS score and nucleosome signal
  metadata$high.tss <- ifelse(metadata$TSS.enrichment > TSS_threshold, 'High TSS score cells', 'Low TSS score cells')
  metadata$nucleosome_group <- ifelse(metadata$nucleosome_signal > nucleosome_threshold, 'nucleosome signal > 4', 'nucleosome signal < 4')
  
  # check for sample_prefixes and add sample names to metadata
  if(!(is.null(sample_names) && is.null(sample_prefixes))) {
    if((length(sample_prefixes) == length(sample_names)) && (length(sample_prefixes) != 0)) {
      for(i in 1:length(sample_prefixes)) {
        if(length(grep(sample_prefixes[i], metadata$cells))==0) {
          stop(paste0('Error: sample prefix ', sample_prefixes[i], ' not found! Please check sample_prefixes!\nExecution halted.'))
        }
        metadata$sample[which(str_detect(metadata$cells, sample_prefixes[i]))] <- sample_names[i]
      }
    } else  if(length(sample_names) == 1){
      metadata$sample <- sample_names
    } else {
      stop('Error: sample_prefixes and sample_names must have the same length!')
    }
  }
    
  
  # add modified meta data to seurat object
  seu@meta.data <- metadata
  
  return(seu)
}

CreateOutputDirectory <- function(save_dir) {
  # creates the default output directory paths in the current working directory using the package 'filenamer'
  # make_path recursively creates directories and subdirectories if not available
  require(filenamer)
  make_path(x = (paste0(save_dir, '/data/')))
  make_path(x = (paste0(save_dir, '/figures/')))
}

Filter_ATAC_QC <- function(seu, cut_PRFs_low = 3000, cut_PRFs_high = 20000, cut_RIPr = 15, cut_BR = 0.05, cut_NS_upper = 4, cut_NS_lower = 0, cut_TSSenr = 2, min_cells_per_feature = 1, cr_version = '2.0') {
  # Filter out low quality reads using selected thresholds - these will change with experiment
  # seu     - seurat object to be filtered
  # cut_XYZ - filter thresholds
  cat('Check for QC metrics...\n')
  # check for presence of QC metrics
  if(sum(c(grepl('pct_reads_in_peaks', colnames(seu@meta.data)),
           grepl('peak_region_fragments', colnames(seu@meta.data)),
           grepl('TSS.enrichment', colnames(seu@meta.data)),
           grepl('blacklist_ratio', colnames(seu@meta.data)),
           grepl('nucleosome_signal', colnames(seu@meta.data)),
           grepl('nucleosome_group', colnames(seu@meta.data)),
           grepl('high.tss', colnames(seu@meta.data)))) != 7) {
    stop('QC values missing. Please first run ComputeQCMetrics on seurat object!')
  }
  
  cat(paste0('Filter data...\n'))
  # Subset seurat object using the defined filter thresholds
  if(compareVersion(cr_version, '2.0') >= 0) {
    seu <- subset(
      x = seu,
      subset = 
        peak_region_fragments > cut_PRFs_low &
        peak_region_fragments < cut_PRFs_high &
        pct_reads_in_peaks > cut_RIPr &
        nucleosome_signal < cut_NS_upper &
        nucleosome_signal > cut_NS_lower &
        TSS.enrichment > cut_TSSenr
    )
  } else if(compareVersion(cr_version, '1.2.0') == 0) {
    seu <- subset(
      x = seu,
      subset = 
        peak_region_fragments > cut_PRFs_low &
        peak_region_fragments < cut_PRFs_high &
        pct_reads_in_peaks > cut_RIPr &
        blacklist_ratio < cut_BR &
        nucleosome_signal < cut_NS_upper &
        nucleosome_signal > cut_NS_lower &
        TSS.enrichment > cut_TSSenr
    )
  } else {
    stop('The given Cell Ranger ATAC version is not supported. Please use Cell Ranger ATAC 1.2.0 or later.')
  }

  
  cat(paste0('Remove all features occurring in less than ', min_cells_per_feature, ' cells...\n'))
  # Extract gene counts
  counts <- GetAssayData(object = seu, slot = "counts")

  # Output a logical vector for every gene on whether the more than zero counts per cell
  nonzero <- counts > 0

  # Sums all TRUE values and returns TRUE if more than min_cells_per_feature TRUE values per gene
  keep_genes <- Matrix::rowSums(nonzero) >= min_cells_per_feature

  # Subset seurat object with keep_genes
  seu <- subset(
    x = seu,
    features = rownames(seu[keep_genes]) 
  )
  return(seu)
  
}

NonLDR_Clustering <- function(seu, reduction = 'lsi', dims = 2:30, save_plot = NULL, plotdir = '', res = 1.0, cluster_only = F, skip_clustering = F, dimplot_group = "seurat_clusters") {
  # Non-linear dimension reduction and clustering
  
  # Performs clustering using LSI components and subsequent tSNE and UMAP embeddings
  # seu           - seurat object
  # reduction     - dimension reduction used for clustering and embeddings
  # dims          - which dimensions to use
  # save_plot     - set to file extension to save plots
  # plotdir       - directory where to save the plots
  # res           - resolution used for clustering
  # cluster_only  - only perform clustering and omit embedding
  
  if(!skip_clustering) {
    cat('Find neighbours and clusters...\n')
    # Find the k nearest neighbors for the dataset
    seu <- FindNeighbors(object = seu, reduction = reduction, dims = dims)
    # Identify cell clusters using a shred nearest neighbor modularity optimization based clustering algorithm
    seu <- FindClusters(object = seu, verbose = FALSE, algorithm = 3, resolution = res)
  }

  if(!is.null(save_plot)) {
    r <- format(res, nsmall = 1)
    cat('Assess clustering performance...\n')
    dist.matrix <- dist(x = Embeddings(object = seu[['lsi']])[, 2:30])
    clusters <- seu$seurat_clusters
    clusters <- as.numeric(as.factor(clusters))
    sil <- cluster::silhouette(clusters, dist.matrix)
    print(summary(sil))
    # cat('Save silhouette plot...\n')
    if(length(clusters) < 2000) {
      h <- 500
    } else {
      h <- 500*length(unique(clusters))/10
    }
    png(file = paste0(plot_dir, 'silhouette_r', r, '.png'), width = 600, height = h)
    plot(sil, col=1:length(unique(clusters)), border=NA, main = paste0('Clustering performance - res ', r))
    dev.off()
  }
  
  if (!cluster_only) {
    # Perform both tSNE and UMAP embeddings
    cat('Run tSNE and UMAP embeddings...\n')
    seu <- RunTSNE(object = seu, reduction = reduction, dims = dims)
    seu <- RunUMAP(object = seu, reduction = reduction, dims = dims)
    
    cat('Create plots...\n')
    #Define arrow positions for UMAP plot (dynamic, no need to change)
    x.u <- min(seu@reductions$umap@cell.embeddings[,1])
    y.u <- min(seu@reductions$umap@cell.embeddings[,2])
    f.u <- 0.9
    # plot UMAP with cell type annotation per cluster
    p_u <- DimPlot(seu, label = T, repel = TRUE, reduction = 'umap', group.by = dimplot_group) +
      ggtitle('') + 
      theme(plot.title = element_blank(), # center title and bold font
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.position = 'none') + # blanking out all axes labels
      geom_segment(aes(x=x.u,y=y.u,xend=x.u,yend=y.u+2*f.u), arrow=arrow(length=unit(0.3,"cm"))) + # introduce sideways arrow
      geom_segment(aes(x=x.u,y=y.u,xend=x.u+2*f.u,yend=y.u), arrow=arrow(length=unit(0.3,"cm"))) + # introduce upward arrow
      annotate("text", x = x.u+f.u, y = y.u-f.u, label = "UMAP 1", size = 3) + # sideways component description
      annotate("text", x = x.u-f.u, y = y.u+f.u, label = "UMAP 2", size = 3, angle = 90) # upward component description
    
    #Define arrow positions for tsne plot (dynamic, no need to change)
    x.t <- min(seu@reductions$tsne@cell.embeddings[,1])
    y.t <- min(seu@reductions$tsne@cell.embeddings[,2])
    f.t <- 3
    # plot UMAP with cell type annotation per cluster
    p_t <- DimPlot(seu, label = T, repel = TRUE, reduction = 'tsne', group.by = dimplot_group) +
      ggtitle('') + 
      theme(plot.title = element_blank(), # center title and bold font
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.position = 'none') + # blanking out all axes labels
      geom_segment(aes(x=x.t,y=y.t,xend=x.t,yend=y.t+2*f.t), arrow=arrow(length=unit(0.3,"cm"))) + # introduce sideways arrow
      geom_segment(aes(x=x.t,y=y.t,xend=x.t+2*f.t,yend=y.t), arrow=arrow(length=unit(0.3,"cm"))) + # introduce upward arrow
      annotate("text", x = x.t+f.t, y = y.t-f.t, label = "tSNE 1", size = 3) + # sideways component description
      annotate("text", x = x.t-f.t, y = y.t+f.t, label = "tSNE 2", size = 3, angle = 90) # upward component description
    
    p <- p_t + p_u
    
    print(p)
    
    if(!is.null(save_plot)) {
      file_types <- c('eps', 'ps', 'tex', 'pdf', 'jpeg', 'tiff', 'png', 'bmp', 'svg')
      if(!(save_plot %in% file_types)) {
        stop('File extension for saving plots not recognized! Please refer to ggsave to select proper extension.')
      }
      
      cat('Save plots...\n')
      ggsave(filename = paste0(plotdir, 'tSNE-UMAP-clusters.', save_plot),
             plot = p,
             width = 10,
             height = 5)
    }
  }


  return(seu)
}

NormLDR <- function(seu, min.cutoff = 'q0', plots = NULL, plotdir = '', dims = 20) {
  # Normalization and linear dimensional reduction
  
  # seu         - seurat object
  # min.cutoff  - either percentile (q5 - 95%) or integer for number of variable features
  # plots       - set to file extension to save component correlation with sequencing depth and elbow plots
  # plotdir     - directory where to save the plots; will save in current working directory if not set
  # dims        - dimensions shown in plots
  
  # Run term frequency inverse document frequency (TF-IDF) normalization
  seu <- RunTFIDF(seu)
  # Find top features based on total number of counts
  seu <- FindTopFeatures(seu, min.cutoff = min.cutoff)
  # Run partial singular value decomposition
  seu <- RunSVD(seu)
  
  # plot correlation coefficients of n dimensions and sequencing depth and save if depth_cor = TRUE
  p <- DepthCor(seu, n = dims)
  q <- ElbowPlot(seu, reduction = 'lsi', ndims = dims)
  plot <- p / q
  print(plot)
  if(!is.null(plots)) {
    ggsave(filename = paste0(plotdir, 'DR-QC.', plots),
           plot = plot,
           width = 6,
           height = 6)
  }
  return(seu)
}

PlotUMAP <- function(seu, cols = NULL, group.by = NULL, f = 0.9, label = FALSE, title = '', coord.tex.size = 3, legend.position = c(.8, .1), fontface_umap_coord = 2, arrow_thickness = 1, reduction = "umap", shuffle = T, ...) {
  x.u <- ceiling(min(seu@reductions$umap@cell.embeddings[,1]))-1
  y.u <- ceiling(min(seu@reductions$umap@cell.embeddings[,2]))-1
  # plot UMAP with coloring per sample
  p <- DimPlot(seu, label = label, repel = TRUE, reduction = reduction, group.by = group.by, cols = cols, shuffle = shuffle, ...) +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'), # center title and bold font
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = legend.position) + # blanking out all axes labels
    geom_segment(aes(x=x.u,y=y.u,xend=x.u,yend=y.u+2*f), arrow=arrow(length=unit(0.3,"cm")), linewidth = arrow_thickness) + # introduce sideways arrow
    geom_segment(aes(x=x.u,y=y.u,xend=x.u+2*f,yend=y.u), arrow=arrow(length=unit(0.3,"cm")), linewidth = arrow_thickness) + # introduce upward arrow
    annotate("text", x = x.u+f, y = y.u-f, label = "UMAP 1", size = coord.tex.size, fontface = fontface_umap_coord) + # sideways component description
    annotate("text", x = x.u-f, y = y.u+f, label = "UMAP 2", size = coord.tex.size, fontface = fontface_umap_coord, angle = 90) # upward component description
  return(p)
}

ReadMultipleSamples.scATAC <- function(data.dir, passed.filters = 500, samples = NULL, fragment_file = '10x', use_peaks = NULL) {
  # read multiple samples using a peak union of all samples
  # data.dir - folder string containing all cellranger output folders for the samples, or string vector with folders for each sample
  # passed.filters - amount of minimum reads per cell that passed filtering
  # samples - vector with sample names
  # use_peaks - Path to bed file with peaks to use. If NULL, peaks from cellranger output files will be merged across all datasets.
  
  # check if data folders have correct formatting and adjust formatting
  if(length(grep('/$', data.dir)) != length(data.dir)) {
    data.dir[which(!grepl('/$', data.dir))] <- paste0(data.dir[which(!grepl('/$', data.dir))], '/')
  }
  if(length(data.dir) == 1) {
    dirs <- list.dirs(path = data.dir, full.names = FALSE, recursive = FALSE)
    dirs <- dirs[(dirs != '') & (dirs != '.snakemake') & (dirs != 'fqs')]
    dirs <- paste0(data.dir, dirs, '/outs/')
  } else {
    dirs <- paste0(data.dir, 'outs/')
  }
  if(is.null(samples)) {
    samples <- basename(dirname(dirs))
  }
  nsample <- length(dirs)
  cat(paste0('Read ', nsample, ' samples...\n'))
  
  # prepare lists to save file reads
  gr <- GRanges()
  meta <- list()
  frags <- list()
  
  # read in data per sample
  for(i in 1:length(dirs)) {
    
    cat(paste0('Sample ', i,': ', samples[i], '\n'))
    cat('Read in peaks and metadata...\n')
    
    # search for input files in sample folder
    files <- list.files(path = dirs[i])
    bed <- grep(pattern = 'peaks\\.bed$', x = files)
    csv <- grep(pattern = 'singlecell\\.csv$', x = files)
    if(fragment_file == '10x') {
      tsv.gz <- grep(pattern = 'fragments\\.tsv\\.gz$', x = files)
    } else {
      files <- c(files, fragment_file)
      tsv.gz <- length(files)
    }
    
    if(is.null(use_peaks) & (length(bed) == 0)) {
      stop(paste0('Error: peaks file missing for this sample!'))
    } else if(length(csv) == 0) {
      stop('Error: metadata file missing for this sample!')
    } else if(length(tsv.gz) == 0) {
      stop('Error: fragment file missing for this sample!')
    }
    
    # read in files
    if(is.null(use_peaks)) {
      peaks <- read.table(
        file = paste0(dirs[i], files[bed]),
        col.names = c("chr", "start", "end")
      )
      gr <- c(gr, makeGRangesFromDataFrame(peaks))
    }

    md <- read.table(
      file = paste0(dirs[i], files[csv]),
      stringsAsFactors = FALSE,
      sep = ",",
      header = TRUE,
      row.names = 1
    )[-1, ] # remove the first row
    # only consider cells called by cellranger
    md <- md[!(md$is__cell_barcode==0), ]
    # perform an initial filtering of low count cells
    md <- md[md$passed_filters > passed.filters, ]
    meta[[i]] <- md
    # create fragment objects
    fr <- CreateFragmentObject(
      path = paste0(dirs[i], files[tsv.gz]),
      cells = rownames(md)
    )
    frags <- c(frags, fr)
  }
  
  # combine all sample peaks by union
  if(is.null(use_peaks)) {
    combined.peaks <- reduce(x = gr)
  } else {
    peaks <- read.table(
      file = use_peaks,
      col.names = c("chr", "start", "end")
    )
    combined.peaks <- makeGRangesFromDataFrame(peaks)
  }
  
  
  # define list for seurat objects
  seus <- list()
  
  # create seurat objects for each sample using the combined peak set
  for(i in 1:length(dirs)) {
    cat(paste0('Sample ', samples[i],':\n'))
    cat('Retrieve counts...\n')
    counts <- FeatureMatrix(
      fragments = frags[[i]],
      features = combined.peaks,
      cells = rownames(meta[[i]])
    )
    cat('Create ChromAssay...\n')
    chrom.assay <- CreateChromatinAssay(counts, fragments = frags[[i]])
    cat('Create Seurat object...\n')
    seu <- CreateSeuratObject(chrom.assay, assay = "peaks", meta.data = meta[[i]])
    seu$sample <- samples[i]
    seus <- c(seus, seu)
  }
  
  # Merge all seurat objects
  cat('Merge samples...\n')
  combined <- merge(
    x = seus[[1]],
    y = seus[c(-1)],
    add.cell.ids = samples
  )
  cat('Data read in successful.\n')
  return(combined)
}
