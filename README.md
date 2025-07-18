# COMICs_Seurat
Tips on how to run Seurat as we do in the lab.

# Set up environment (if running on cluster)

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

plan("multicore", workers = 2)  # Adjust cores needed - With large objects less is more!

options(future.globals.maxSize = 256 * 1024^3) # Correct RAM - In this case, more is more!

setwd("~/My/directory") # Set working directory

````

## Read, QC and process

````

# If you already have a pre-processed object in rds

seurat <- readRDS("my_object.rds")

# If you have 10X Output

# Load data - Do like this for unfiltered data from Cell Ranger ####

raw_counts <- Read10X(data.dir = "10X_Out_Dir/")

# Remove empty droplets
e.out <- emptyDrops(raw_counts, lower = 200)

# Retain only droplets with FDR < 0.01 (or < 0.05 if you want to be more inclusive)
e.out$FDR[is.na(e.out$FDR)] = 1
valid_barcodes <- which(e.out$FDR < 0.01)

filtered_counts <- raw_counts[, valid_barcodes]

# Create Seurat Object with Valid Barcodes and filtered
seurat<- CreateSeuratObject(counts = filtered_counts, min.cells = 3, min.features = 200)

#### Initial QC ####

# Check mitochondrial gene expression
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter based on the plots -> Subjective filters, consider the characteristics of your data!
# Some guide lines: [200,500] < nFeature_RNA < [2500,5000] ; percent.mt < 5 for normal cells and < 20 for tumor cells

seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#### Normalization and Scaling ####

# For this step we will use the SCTransform function from seurat. But first, we will run the standard pipeline
# and check if there is separation in the PCA due to a specific variable. If so, we also select that var to be
# regressed in the SCTransform

#Normalize
seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
seurat <- FindVariableFeatures(object = seurat, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data - Takes a while!
seurat <- ScaleData(object = seurat, features = VariableFeatures(seurat))

#Run PCA and check for possible separations based on specific variables

seurat <- RunPCA(object = seurat, features = VariableFeatures(object = seurat)) 
DimPlot(seurat, reduction = "pca", group.by = "percent.mt") + NoLegend()

# Run SCTransform. Adjust vars to regress as needed

seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)

#Run PCA and Determine Dimensions for 90% Variance. Correct Batch effect with Harmony if needed

seurat <- RunPCA(object = seurat, features = VariableFeatures(object = seurat))
seurat <- RunHarmony(seurat, group.by.vars = c("Variables_of_Interest_For_Batch_Effect"), plot_convergence = F) # To check if you have batch effect you may need to do the reach the UMAP Step and plot according to varaibles like sample origin

stdev <- seurat@reductions$harmony@stdev #Change to pca if you didn't run harmony!
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

#Find Neighbors + Find CLusters 

seurat <- FindNeighbors(seurat, dims = 1:PCNum, verbose = TRUE, reduction = "harmony")
seurat <- FindClusters(seurat, verbose = TRUE, resolution = 0.4)
seurat <- RunUMAP(seurat, dims = 1:PCNum, verbose = TRUE, reduction = "harmony", n.components = 2L)

# Save processed object
saveRDS(seurat,"my_object_processed.rds")

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

## Optional QC - Doublet Finder

DoubletFinder generates artifical doublets and then does a PCA to try to identify cells that are close to this artificial doublets. To run DoubletFinder you need to run everything until the FindClusters step. Adjust the expected doubled ratio based on this thread https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76#issuecomment-624825029

````
library(DoubletFinder)

# Estimate expected doublet rate 
nExp <- round(0.050 * ncol(seurat))

# Run DoubletFinder
seurat <- doubletFinder(seurat,
                              PCs = 1:PCNum,
                              pN = 0.25,
                              pK = 0.005,  # you cab try paramSweep to optimize this
                              nExp = nExp,
                              reuse.pANN = NULL,
                              sct = TRUE)
  
# Remove predicted doublets
seurat <- subset(seurat, subset = DF.classifications_0.25_0.005_559 == "Singlet") # Note: DF will generate a column with a different name in your metadata. Check the name!
  
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



## Adjusting resolution 

Changing the number of clusters can be done by adjusting the resolution. Clustree allows to have an idea of how the clusters change with different resolutions.

````
library(clustree)

seurat <- FindClusters(seurat, verbose = TRUE, resolution = c(0, 0.2, 0.4, 0.6, 0.8, 1))

clustree(seurat)
````

## Subseting

Sometimes you will want to consider only a subset of cells. This can be easily done:

````
fibro_subset <- subset(seurat, subset = celltype == "Fibroblast")
````

Don't forget to re-do the normalization steps on the new subset!

## Sample Integration

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
