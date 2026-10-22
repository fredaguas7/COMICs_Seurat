# COMICs_Seurat
Tips on how to run Seurat as we do in the lab.

# Set up environment

You can load my environment directly (as of 22 October 2026) from the yml file in this repo

````
conda env create -f seurat_FA.yml
````

You can create your own and install the packages like this (not complete, install what you need)

```
# Create env 
conda create -n seurat 

# Install r and packages 

conda install -c conda-forge r-base=4.4.2 r-essentials

# Install the packages also via conda to avoid clashes 

conda install -c bioconda -c conda-forge r-seurat :bioconductor-glmgampoi bioconductor-dropletutils 
conda install -c conda-forge r-patchwork r-dplyr r-harmony r-clustree
conda install genomedk::r-doubletfinder

conda install conda-forge::r-ggplot2
conda install conda-forge::r-future
```
# Running Seurat

## Libraries and set-up

```
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(ggplot2)
library(future)
library(glmGamPoi)
library(DropletUtils)
library(Matrix)
library(SoupX)
library(DoubletFinder)
library(clustree)
library(scater) # For MAD
library(SingleCellExperiment) # For MAD
library(RColorBrewer)
library(edgeR)
library(data.table)
library(rstatix) #Stats
library(coin) #Plots
library(ggpubr) # Plots
library(scales) # Plots 
library(ggalluvial) # Sankey
library(regioneR) # Permutations
library(tidyr)
library(ggridges)
library(AUCell)
library(pheatmap) # Heatmap

plan("multicore", workers = 2)  # Adjust cores needed - With large objects less is more!

options(future.globals.maxSize = 256 * 1024^3) # Correct RAM - In this case, more is more!

setwd("~/My/directory") # Set working directory

````

## Read, SoupX, EmptyDrops and DoubletFinder

This is being relatively strict, remember to adjust the parameters if needed. Basically, read the Star_Solo output (raw), adjust ambient RNA with SoupX, remove empty drops and doublets.

````

base_dir <- "/My/Output/From/STAR_Solo"

out_dirs <- list.dirs(
  path = base_dir,
  recursive = FALSE,
  full.names = TRUE
)

out_dirs <- out_dirs[grepl("\\.out$", out_dirs)] # Assuming base StarSolo names

#Load raw for SoupX and emptyDrops! 
gene_dirs <- file.path(out_dirs, "Gene", "raw")

# Read and creat seurat object

star_counts <- lapply(gene_dirs, function(d) {
  ReadSTARsolo(data.dir = d)
})

# Correct sample names, assuming base StarSolo names
sample_names <- sub("\\_Solo.out$", "", basename(dirname(dirname(gene_dirs))))
names(star_counts) <- sample_names

# Remove MT/RB genes for emptyDrops as recommended

ref_genes <- rownames(star_counts[[1]])
mt_genes  <- grep("^MT-", ref_genes, value = TRUE, ignore.case = TRUE)
rb_genes  <- grep("^RP[LS]", ref_genes, value = TRUE, ignore.case = TRUE)
qc_exclude <- c(mt_genes, rb_genes)


# Remove empty droplets with emptyDrops
filtered_counts <- lapply(names(star_counts), function(sample) {
  
  message("Running Strict emptyDrops (MT/Ribo excluded) for ", sample)
  
  # emptyDrops requires integer counts
  counts <- round(star_counts[[sample]])
  
  # Identify indices of genes to keep for the significance calculation
  keep_genes <- !rownames(counts) %in% qc_exclude
  
  # Run emptyDrops ONLY on the non-MT/non-Ribo subset
  e.out <- DropletUtils::emptyDrops(counts[keep_genes, ], lower = 200)
  
  # Correct NAs to empty droplets
  e.out$FDR[is.na(e.out$FDR)] <- 1
  
  # Use the Strict FDR threshold - 0.001
  valid_barcodes <- which(e.out$FDR < 0.001)
  
  message(sample, ": Retained ", length(valid_barcodes), " high-confidence droplets.")
  
  # Return the FULL matrix (all genes) for the valid barcodes
  counts[, valid_barcodes, drop = FALSE]
})

names(filtered_counts) <- names(star_counts)

# ── MAD-based QC on emptyDrops output ────────────────────────────────────────

filtered_counts <- lapply(names(filtered_counts), function(sample) {
  
  message("Running MAD QC for ", sample)
  
  mat     <- filtered_counts[[sample]]
  sce     <- SingleCellExperiment(assays = list(counts = mat))
  is.mito <- grepl("^MT-", rownames(sce))
  qc      <- perCellQCMetrics(sce, subsets = list(mito = is.mito))
  
  low_lib   <- isOutlier(qc$sum,                log = TRUE,  type = "lower",  nmads = 3)
  low_gene  <- isOutlier(qc$detected,           log = TRUE,  type = "lower",  nmads = 3)
  high_mito <- isOutlier(qc$subsets_mito_percent, log = FALSE, type = "higher", nmads = 3) |
    qc$subsets_mito_percent > 25
  
  discard <- low_lib | low_gene | high_mito
  
  message(sample, ": low_lib=",        sum(low_lib),
          " | low_gene=",      sum(low_gene),
          " | high_mito=",     sum(high_mito),
          " | total_removed=", sum(discard),
          " | keeping=",       sum(!discard))
  
  mat[, !discard, drop = FALSE]
})

names(filtered_counts) <- names(star_counts)

# Number of cells per sample for doublet rate
cell_counts_df <- data.frame(
  sample = names(filtered_counts),
  cell_count = sapply(filtered_counts, ncol)
)

print(cell_counts_df)

seurat_list <- lapply(names(filtered_counts), function(sample) {
  message("Processing: ", sample)
  
  # --- 1. SOUPX ---
  raw_mat <- star_counts[[sample]]
  filt_mat <- filtered_counts[[sample]]
  sc <- SoupChannel(raw_mat, filt_mat) # Takes the raw vs the empty droplets to build the soup
  
  # Use a temporary object for SoupX clustering
  tmp <- CreateSeuratObject(counts = filt_mat)
  tmp <- NormalizeData(tmp, verbose = FALSE)
  tmp <- FindVariableFeatures(tmp, verbose = FALSE)
  tmp <- ScaleData(tmp, verbose = FALSE)
  tmp <- RunPCA(tmp, verbose = FALSE)
  tmp <- FindNeighbors(tmp, dims = 1:20, verbose = FALSE)
  tmp <- FindClusters(tmp, resolution = 0.5, verbose = FALSE, algorithm=4)
  tmp <- RunUMAP(tmp, dims = 1:20, verbose = FALSE)
  
  # Ensure clusters are a named vector, not a DF
  sc_clusters <- setNames(as.character(tmp$seurat_clusters), colnames(tmp))
  sc <- setClusters(sc, sc_clusters)
  sc <- autoEstCont(sc, doPlot = FALSE) # Estimate contmaination
  out <- adjustCounts(sc, roundToInt = TRUE) # Clean up
  
  # --- 2. CREATE CLEAN OBJECT ---
  obj <- CreateSeuratObject(counts = out, min.cells = 3, min.features = 200, project = sample)
  obj$sample <- sample
  obj$soup_rho <- sc$fit$rhoEst # Save the contamination proportion
  
  # --- 3. PREP FOR DOUBLETFINDER ---
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE, algorithm=4)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  
  # ---  Find Optimal pK ---
  sweep.res.list <- paramSweep(obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  bcmvn$pK <- as.numeric(as.character(unlist(bcmvn$pK)))
  
  optimal.pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  
  # Take clusters to estimate the homotypic (same celltype in a doublet) proportion
  annotations <- obj@meta.data$seurat_clusters

  annotations <- as.vector(annotations)
  annotations <- as.character(annotations)
  homotypic.prop <- modelHomotypic(annotations)           
  
  # May need to adjust this rate
  # Formula: (cells / 1000) * 0.008
  n_cells <- nrow(obj@meta.data)
  dynamic_rate <- (n_cells / 1000) * 0.008
  
  # Calculate expected number of doublets
  nExp_poi <- round(dynamic_rate * n_cells)
  nExp_poi.adj <- max(round(nExp_poi * (1 - homotypic.prop)), 1)
  
  
  # Run DoubletFinder
  obj <- doubletFinder(obj, PCs = 1:20, pN = 0.25, pK = optimal.pk, 
                       nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
  
  # Robustly pull the classification column
  df_col_name <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  obj$doublet_status <- as.character(obj@meta.data[[df_col_name]])
  
  # Rename cells ONLY ONCE at the end
  obj <- RenameCells(obj, add.cell.id = sample)
  
  return(obj)
})
  

names(seurat_list) <- names(star_counts)

#Sanity checks
lapply(seurat_list, dim)
lapply(seurat_list, function(x) head(colnames(x)))
lapply(seurat_list, function(x) head(rownames(x)))


# Merge into one object
seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  project = "RNA"
)

seurat <- JoinLayers(seurat)

# Filter doublets out after merging
seurat <- subset(seurat, subset = doublet_status == "Singlet")

# Filter genes not detected in at least 20 cells across the full dataset
counts_mat   <- GetAssayData(seurat, layer = "counts")
genes_keep   <- rowSums(counts_mat > 0) >= 20
seurat       <- seurat[genes_keep, ]

````

## In QC remember to adjust thresholds based on your samples!

Regarding %mt, for cancer I would not go lower than 20/25%. For normal cells 5% will be enough.

````

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, "^MT")

# QC metrics (violin plot)

p <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =1, pt.size = 0)
ggsave("QC_VlnPlot_Sample.pdf", plot = p , width = 15, height = 8)

p <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "condition", ncol =3, pt.size = 0)
ggsave("QC_VlnPlot_Condition.pdf", plot = p , width = 15, height = 8)

# filter

seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25 ) # Adjust as needed

````

This QC scatterplot looks at complexity by sample


````

QC_scatterplot_by_sample <- function(seurat_obj, x_var = "nCount_RNA", y_var = "nFeature_RNA", col_var = "percent.mt") {
  
  p <- ggplot(seurat_obj@meta.data, aes(.data[[x_var]], y = .data[[y_var]], color = .data[[col_var]])) +
    geom_point(alpha = 0.5, size = 0.05) +
    theme_bw() +
    facet_wrap(~sample, ncol = 5) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "UMI Count (log10)", 
      y = "Gene Count (log10)", 
      color = col_var
    ) +
    theme(
      axis.title = element_text(color = "black", size = 10),
      axis.text = element_text(colour = "black", size = 8),
      axis.text.x = element_text(angle = -45, hjust = 0),
      text = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 8),
      strip.background = element_rect(fill = "white", color = "black"),
      legend.position = "top",
      legend.key.width = unit(1.2, "cm")
    )
  
  # Apply color scale based on col_var
  if (col_var == "celltype") {
    p <- p + scale_color_manual(values = celltype_colors)
  } else {
    p <- p + scale_colour_gradientn(
      colours = c("#2ECC40", "#FFFF00", "#FF4136")
    )
  }
  
  return(p)
}


# Adjust sizes as needed
p <- QC_scatterplot_by_sample(seurat)
ggsave("Complexity_By_Sample.pdf", p, width = 15, height = 20)

````
## For normalization we stick with SCT without layer splitting

Do run your "normal" seurat pipeline first to check the PCA for %mt driving separation (almost always the case, if so regress for %mt in SCT)

````


#Normalize
seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
seurat <- FindVariableFeatures(object = seurat, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
seurat <- ScaleData(object = seurat, features = VariableFeatures(seurat))

#Run PCA and check for possible separations based on specific variables (percent mt)

seurat <- RunPCA(object = seurat, features = VariableFeatures(object = seurat))


p <- DimPlot(seurat, reduction = "pca", group.by = "percent.mt") + NoLegend()
ggsave("PCA_percentmt.pdf", plot = p , width = 10, height = 8)

p <- DimPlot(seurat, reduction = "pca", group.by = "nCount_RNA") + NoLegend()
ggsave("PCA_Count.pdf", plot = p , width = 10, height = 8)


#### RUN SCT ####

seurat <- SCTransform(seurat, vars.to.regress = c("nCount_RNA","percent.mt"), verbose = TRUE, conserve.memory = TRUE)

#Run PCA 

seurat <- RunPCA(object = seurat, features = VariableFeatures(object = seurat))

````

## Batch effect correction with harmony

I usually take enough PCs to explain 90% of variance. Might be worth checking the elbow plot anyway 

````

seurat <- RunHarmony(seurat, group.by.vars = c("sample"), plot_convergence = F) 

stdev <- seurat@reductions$harmony@stdev 
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}

#Confirm PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)

PCNum

p <- ElbowPlot(seurat,ndims=50)
ggsave("ElbowPlot.pdf", plot = p , width = 10, height = 8)
````

## Clustering

Currently running Leiden as it seems to me recommended over Louvain. Run multiple resolution and check clustree. Adjust based on your data

````

seurat <- FindNeighbors(seurat, dims = 1:PCNum, verbose = TRUE, reduction = "harmony")

# Algorithm = 4 -> Leiden

seurat <- FindClusters(seurat, algorithm = 4, verbose = TRUE, resolution = c(0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5))

p <- clustree(seurat)
ggsave("Clustree.pdf", plot = p , width = 10, height = 8)


#### UMAPs ####

# Run UMAP
seurat <- RunUMAP(seurat, dims = 1:PCNum, verbose = TRUE, reduction = "harmony", n.components = 2L)

# Set the identity based on the best res
Idents(seurat) <- "SCT_snn_res.0.2"

````

## Visualize 

````
# UMAP for the seurat clusters

p <- DimPlot(object = seurat_filtered, reduction = "umap", label = TRUE, pt.size = 0.5, group.by="seurat_clusters") + NoLegend()

ggsave("UMAP.pdf", plot = p , width = 10, height = 8)

# UMAP for other variable, like cell type

p <- DimPlot(object = seurat_filtered, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "celltype") + NoLegend()

ggsave("UMAP_celltype.pdf", plot = p , width = 10, height = 8)

# Dotplots to visualize expression of some genes. You can also use VlnPlot
# Other nice tools in this plot include group.by and split.by

features_to_plot <- c('Gene_A','Gene_B','Gene_C')

#For the DotPlot is recommended to visualize the original counts! assay='RNA'
p <- DotPlot(object = seurat, features = features_to_plot, cols = "RdBu", 
             group.by = "Clusters", assay='RNA') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("Dotplot.pdf"), plot = combined_plot, width = 10, height = 8)


````

## Find markers 

````
# Find Markers Based on clusters

# Do you want only positive markers? Only LogFC > 1?
# How strict do you want to be regarding the % of cells expressing the marker?

seurat <- PrepSCTFindMarkers(seurat)

seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, assay = "SCT",
                                          logfc.threshold = 1,min.pct=0.2)

# Eles têm esta info na vignette do SCT mas não devemos sempre re-normalizar quando se faz um subset? Acho que não faz sentido ter isto
# If running on a subset of the original object after running PrepSCTFindMarkers(), FindMarkers() should be invoked with recorrect_umi = FALSE to use the existing corrected counts

write.table(seurat.markers, "All_markers.csv")

# Find Markers based on other variable - ex, celltype

Idents(seurat) <- seurat_filtered@meta.data$celltype

seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE,
                                          logfc.threshold = 1,min.pct=0.2) 

write.table(seurat.markers, "All_markers_celltype.csv")
````

## Assign your cluster names based on markers you saw or a dotplot of known markers. Example:

````

clusters <- FetchData(seurat, vars = "SCT_snn_res.0.2")

cell_map <- c(
  "1" = "Ductal 1",
  "2" = "T",
  "3" = "Ductal 2",
  "4" = "Endothelial",
  "5" = "Stellate",
  "6" = "Macrophage",
  "7" = "B",
  "8" = "Fibroblasts",
  "9" = "Acinar",
  "10" = "Ductal 1",
  "11" = "Endocrine",
  "12" = "Ductal 1",
  "13" = "Mast",
  "14" = "T"
)

celltypes <- cell_map[as.character(clusters$SCT_snn_res.0.2)] 

names(celltypes) <- rownames(clusters)

seurat[["celltype"]] <- celltypes

Idents(seurat) <- "celltype"
````

## Check if your cell type is actually the main driver of clustering

Plot the importance of certain columns to the PCA to check that the driver is the cell type: 

````

plot_pc_variance <- function(seurat_obj, variables, num_pcs = 30, reduction = "pca") {
  
  # 1. Check if variables exist in metadata
  missing <- variables[!variables %in% colnames(seurat_obj@meta.data)]
  if(length(missing) > 0) {
    stop(paste("The following variables are not in the metadata:", paste(missing, collapse = ", ")))
  }
  
  # 2. Extract PCA embeddings
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in Seurat object."))
  }
  
  pcs <- Embeddings(seurat_obj, reduction = reduction)[, 1:num_pcs]
  
  # To avoid reductions with different sizes -> TEs
  cells_in_reduction <- rownames(pcs)
  
  # 3. Calculate R-squared for each variable
  pcr <- lapply(variables, function(x) {
    # Extract metadata column
    meta_val <- seurat_obj@meta.data[cells_in_reduction, x, drop = TRUE]
    
    # Run linear model across all selected PCs
    # lm(Matrix ~ vector) fits a separate model for each column of the matrix
    fit <- summary(lm(pcs ~ meta_val))
    
    # Extract Adjusted R-squared
    r2 <- sapply(fit, function(f) f$adj.r.squared)
    
    data.frame(
      Variable = x, 
      PC = seq_along(r2), 
      R2 = r2
    )
  }) %>% bind_rows()
  
  # 4. Generate Plot
  ggplot(pcr, aes(x = PC, y = R2, color = Variable)) +
    geom_line(linewidth = 1) + 
    geom_point(size = 2) +
    theme_bw() +
    scale_x_continuous(breaks = seq(1, num_pcs, by = 2)) +
    labs(
      title = "Principal Component Variance Analysis",
      subtitle = paste("Reduction:", reduction),
      x = "Principal Component",
      y = "Adjusted R-squared (Variance Explained)",
      color = "Metadata"
    ) +
    theme(
      text = element_text(face = "bold"),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}


target_vars <- c("nCount_RNA", "percent.mt", "sample", "condition", "celltype")
p <- plot_pc_variance(seurat, variables = target_vars)
ggsave("PC_Variance_Explanation.pdf", p, width = 10, height = 6)

p <- plot_pc_variance(seurat, variables = target_vars, reduction = "harmony")
ggsave("PC_Variance_Explanation_Harmony.pdf", p, width = 10, height = 6)

````


## Subseting

Sometimes you will want to consider only a subset of cells. This can be easily done:

````
fibro_subset <- subset(seurat, subset = celltype == "Fibroblast")
````

Don't forget to re-do the normalization steps on the new subset!

## Pseudobulk

Due to previous problems with imbalances in cell numbers we are now running permutation for pseudobulk.

````

run_pseudobulk_permutations <- function(
    seurat_obj,
    assay = "TE",
    layer = "counts",
    celltype_col = "celltype",
    sample_col = "sample",
    condition_col = "condition",
    contrast_levels = c("Normal", "Tumor"), # c(Reference, Target)
    cell_types = NULL,
    n_permutations = 1000,
    stability_threshold = 0.90,
    downsample_pct = 0.90,
    min_cells_per_sample = 10,
    min_samples_per_group = 3,
    output_dir = NULL,
    seed = 42,
    max_fdr = 0.05,
    min_logfc = 1,
    save_all_permutations = TRUE
) {
  set.seed(seed)
  
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Extract count matrix
  if (packageVersion("Seurat") >= "5.0.0") {
    raw_counts <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  } else {
    raw_counts <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer)
  }
  
  cell_meta <- seurat_obj@meta.data
  
  if (is.null(cell_types)) {
    cell_types <- unique(cell_meta[[celltype_col]])
    cell_types <- cell_types[!is.na(cell_types)]
  }
  
  results_list <- list()
  
  for (cell_type in cell_types) {
    safe_cell_name <- gsub("[/ ]", "_", cell_type)
    message(paste("\n--- Running Permutations for Cell Type:", cell_type, "---"))
    
    sub_meta <- cell_meta[cell_meta[[celltype_col]] == cell_type & !is.na(cell_meta[[celltype_col]]), ]
    
    valid_barcodes <- intersect(rownames(sub_meta), colnames(raw_counts))
    if (length(valid_barcodes) == 0) {
      message(paste("Skipping", cell_type, "- no matching barcodes between metadata and count matrix."))
      next
    }
    sub_meta <- sub_meta[valid_barcodes, , drop = FALSE]
    
    sample_counts <- table(sub_meta[[sample_col]])
    valid_samples <- names(sample_counts[sample_counts >= min_cells_per_sample])
    
    sample_cond_df <- unique(sub_meta[sub_meta[[sample_col]] %in% valid_samples, c(sample_col, condition_col)])
    cond_vec <- setNames(sample_cond_df[[condition_col]], sample_cond_df[[sample_col]])[valid_samples]
    conditions <- factor(cond_vec, levels = contrast_levels)
    
    samples_per_group <- table(conditions)
    if (any(samples_per_group < min_samples_per_group)) {
      counts_str <- paste(names(samples_per_group), samples_per_group, sep = "=", collapse = ", ")
      message(paste("Skipping", cell_type, "- requires >=", min_samples_per_group,
                    "valid samples per group. Found:", counts_str))
      next
    }
    
    abs_min <- min(sample_counts[valid_samples])
    target_n_cells <- max(1, floor(abs_min * downsample_pct))
    
    cells_by_sample <- split(rownames(sub_meta), sub_meta[[sample_col]])[valid_samples]
    all_features  <- rownames(raw_counts)
    
    iter_names    <- paste0("iter_", seq_len(n_permutations))
    sig_tracker   <- matrix(0, nrow = length(all_features), ncol = n_permutations, dimnames = list(all_features, iter_names))
    logfc_tracker <- matrix(NA_real_, nrow = length(all_features), ncol = n_permutations, dimnames = list(all_features, iter_names))
    fdr_tracker   <- matrix(NA_real_, nrow = length(all_features), ncol = n_permutations, dimnames = list(all_features, iter_names))
    
    # ------------------ Permutation Loop ------------------
    for (iter in 1:n_permutations) {
      pb_list <- lapply(valid_samples, function(s) {
        sampled_barcodes <- sample(cells_by_sample[[s]], size = target_n_cells, replace = FALSE)
        Matrix::rowSums(raw_counts[, sampled_barcodes, drop = FALSE])
      })
      
      pb_counts <- do.call(cbind, pb_list)
      colnames(pb_counts) <- valid_samples
      
      # edgeR Pipeline
      y <- edgeR::DGEList(counts = pb_counts, group = conditions)
      
      keep <- edgeR::filterByExpr(y, group = conditions)
      y <- y[keep, , keep.lib.sizes = FALSE]
      
      y <- edgeR::calcNormFactors(y)
      
      if (nrow(y) < 10) next
      
      design_matrix <- stats::model.matrix(~ 0 + conditions)
      colnames(design_matrix) <- levels(conditions)
      
      contrast_formula <- paste0(levels(conditions)[2], " - ", levels(conditions)[1])
      contrast_matrix  <- limma::makeContrasts(contrasts = contrast_formula, levels = design_matrix)
      
      y <- edgeR::estimateDisp(y, design_matrix)
      fit <- edgeR::glmQLFit(y, design_matrix)
      qlf <- edgeR::glmQLFTest(fit, contrast = contrast_matrix)
      
      res <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
      
      tested_genes <- rownames(res)
      logfc_tracker[tested_genes, iter] <- res$logFC
      fdr_tracker[tested_genes, iter]   <- res$FDR
      
      sig_genes <- tested_genes[res$FDR < max_fdr & abs(res$logFC) >= min_logfc]
      sig_tracker[sig_genes, iter] <- 1
    }
    
    # ------------------ Results Aggregation ------------------
    selection_freq <- rowMeans(sig_tracker)
    median_logFC   <- apply(logfc_tracker, 1, stats::median, na.rm = TRUE)
    median_FDR     <- apply(fdr_tracker, 1, stats::median, na.rm = TRUE)
    
    # Directional consistency counts
    n_up       <- rowSums(logfc_tracker > 0, na.rm = TRUE)
    n_down     <- rowSums(logfc_tracker < 0, na.rm = TRUE)
    n_tested   <- rowSums(!is.na(logfc_tracker))
    
    concordance <- ifelse(n_tested > 0, pmax(n_up, n_down) / n_tested, NA_real_)
    direction   <- ifelse(n_up == n_tested & n_tested > 0, "Up",
                          ifelse(n_down == n_tested & n_tested > 0, "Down", "Inconsistent"))
    
    permutation_summary <- data.frame(
      Feature               = all_features,
      Selection_Frequency   = selection_freq,
      Median_logFC          = median_logFC,
      Median_FDR            = median_FDR,
      N_Up                  = n_up,
      N_Down                = n_down,
      Direction_Concordance = round(concordance, 4),
      Direction             = direction,
      stringsAsFactors      = FALSE
    )
    
    permutation_summary <- permutation_summary[order(permutation_summary$Selection_Frequency, decreasing = TRUE), ]
    
    # Subset of robust, directionally consistent DEGs for reference
    stable_degs <- permutation_summary[
      permutation_summary$Selection_Frequency >= stability_threshold &
        permutation_summary$Direction != "Inconsistent",
    ]
    
    cell_type_results <- list(
      summary     = permutation_summary,
      stable_degs = stable_degs
    )
    
    if (save_all_permutations) {
      cell_type_results$all_permutations <- list(
        logFC = logfc_tracker,
        FDR   = fdr_tracker
      )
    }
    
    results_list[[cell_type]] <- cell_type_results
    
    if (!is.null(output_dir)) {
      csv_path <- file.path(output_dir, paste0(safe_cell_name, "_Permutation_Results.csv"))
      utils::write.csv(permutation_summary, file = csv_path, row.names = FALSE)
      
      rds_path <- file.path(output_dir, paste0(safe_cell_name, "_Permutation_Results.rds"))
      saveRDS(cell_type_results, file = rds_path)
    }
  }
  
  return(results_list)
}



# Run the permutation analysis
perm_results <- run_pseudobulk_permutations(
  seurat_obj           = seurat,
  assay                = "RNA",
  celltype_col         = "celltype",
  sample_col           = "sample",
  condition_col        = "condition",
  contrast_levels      = c("Normal", "Tumor"),
  cell_types = "Ductal 2",
  n_permutations       = 1000,
  stability_threshold  = 0.90,
  downsample_pct       = 0.90,
  min_cells_per_sample = 10,
  min_samples_per_group= 3,
  output_dir           = "Test_Pseudobulk_Permutations"
)
````

## Sample Integration [WORK IN PROGRESS]

For sample integration you should process each object individually (until RunUMAP) and then integrate

````
# Don't forget to keep a column in the metadata to know what is what

seurat_ctrl$orig.ident <- "control"
seurat_condition$orig.ident <- "condition"

data_list <- list(ctrl = seurat_ctrl, cond = seurat_condition)
features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000)
data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT",
                                  anchor.features = features)
seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT") # Default assay will now be 'integrated'. When FindingMarkers don't forget to use assay="SCT" or change the DefaultAssay to SCT!
  
````
