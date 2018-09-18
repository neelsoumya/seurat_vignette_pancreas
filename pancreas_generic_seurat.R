######################################################################
# Sample workflow in Seurat generic functions
#
# Usage: nohup R --no-save < pancreas_generic_seurat.R
# OR
# R # then
# source("pancreas_generic_seurat.R")
#
# Installation:
#   install.packages('Seurat')
#
# Adapted from
# https://www.dropbox.com/s/irx9qriblxqar14/pancreas.R?dl=1
#
###############################################################


###########################
# Load libraries
###########################

library(Seurat)
library(Matrix)

######################################################################
# Sample workflow for leading in files from GEO
######################################################################
# data downloaded from http://satijalab.org/seurat/get_started.html
# and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133


###########################
# Load data
###########################

str_dir_files = "~/Downloads/GSE84133_RAW/"
str_output_Robject = "pancreas.Robj"

all.files <- list.files(str_dir_files)
all.data <- data.frame()
for (i in all.files[1:4]) {
  dat <- read.csv(paste(str_dir_files, i, sep = ""))
  all.data <- rbind(all.data,data.frame(dat))
  print(i)
}


###########################
# Data transformations
###########################

new.data <- t(all.data[, c(-1, -2, -3)]) # transposing
colnames(new.data) <- all.data[, 1]
pancreas.data <- new.data
pancreas.md <- all.data[, 2:3]
rownames(pancreas.md) <- all.data[, 1]

######################################################
# Quality control and preprocessing
# 1) Filter cells to get a number of cells after QC
#       a. cells expressing fewer than 500 genes per cell
#       b. genes expressed in less than 3 cells were removed from the analysis
#
# 2) Normalize
# 3) Find variable genes
# 4) Remove unwanted variation (use a linear model to remove unwanted variation)
# 5) Run PCA on variable genes and tSNE also 
# 6) Cluster the data on these principal components
# 7) Find markers of these clusters
######################################################

pancreas.data <- Matrix(pancreas.data, sparse = T)
pancreas <- CreateSeuratObject(raw.data = pancreas.data, min.cells = 3)
pancreas <- AddMetaData(pancreas, metadata = pancreas.md)

# Filter cells to get a number of cells after QC
#       a. cells expressing fewer than 500 genes per cell
#       b. genes expressed in less than 3 cells were removed from the analysis
pancreas <- FilterCells(pancreas, subset.names = "nGene", low.thresholds = 500,  high.thresholds = Inf)
pancreas <- NormalizeData(pancreas)

# Find variable genes
pancreas <- FindVariableGenes(pancreas, x.low.cutoff = 0.1)

# Remove unwanted variation (use a linear model to remove unwanted variation)
pancreas <- ScaleData(pancreas, vars.to.regress = c("orig.ident", "nUMI"), genes.use = pancreas@var.genes, model.use = "negbinom")

# Run PCA on variable genes and tSNE also 
pancreas <- RunPCA(pancreas, pcs.compute = 30, weight.by.var = FALSE)

# Get scree plot
PCElbowPlot(object = pancreas)

# Get jackStraw plots
pancreas <- JackStraw(object = pancreas, num.replicate = 100, display.progress = FALSE)

pancreas <- RunTSNE(pancreas, dims.use = 1:19, do.fast = T)

# Cluster the data on these principal components
pancreas <- FindClusters(pancreas, reduction.type = "pca", dims.use = 1:19, save.SNN = T)

# color by cluster ID, annotated cluster from the manuscript, or batch
# Can switch the identity class using SetAllIdent if desired

###########################
# tSNE plots
###########################

TSNEPlot(pancreas, do.label = TRUE)
TSNEPlot(pancreas, group.by = "assigned_cluster")
TSNEPlot(pancreas, group.by = "orig.ident")


###########################
# Find marker genes
###########################

# Find Markers of ductal cell subcluster, using the negative binomial test
# only test genes with a 20% difference in detection rate to speed-up (optional)
# Find markers of these clusters
ductal.markers <- FindMarkers(pancreas, ident.1 = 5, ident.2 = 12, test.use = "negbinom", min.diff.pct = 0.2)


###########################
# Visualize markers
###########################

# Visualize canonical and new markers
FeaturePlot(pancreas, c("GCG", "INS","TFF1","PRSS1","VGF","TRIB3","DDR1","CRYBA2","SLC30A8"),
            cols.use = c("lightgrey","blue"), nCol = 3)


###########################
# Save output
###########################

# Can save the object for future loading
save(pancreas, file = str_output_Robject)
