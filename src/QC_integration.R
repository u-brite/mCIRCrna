#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(xlsx)
library(ggplot2)

#For 5 m/o mice
#Load the 5m/o filtered matrices from MyoAtlas paper
months5_filtered <- Read10X_h5(filename = "5mo_filtered_feature_bc_matrix.h5", use.names = TRUE)

#Creating the R seurat object for 5 month old mice
fivemonth <- CreateSeuratObject(counts = months5_filtered, project = "A", min.cells = 3, min.features = 200)

#Checking percentage of reads that map to mitochondrial genome
fivemonth[["percent.mt"]] <- PercentageFeatureSet(fivemonth, pattern = "^MT-")

#Vizualizing the QC metrics in Violin plot
VlnPlot(fivemonth, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#Visualizing feature-feature relationship using feature scatterplot
featurescatterplot <- FeatureScatter(fivemonth, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot(featurescatterplot)

#Filtering the reads with features < 200 and >3200
fivemonth <- subset(fivemonth, subset = nFeature_RNA > 200 & nFeature_RNA < 3200)

#Normalizing the Data 
fivemonth <- NormalizeData(fivemonth)

#identifying highly variable features
fivemonth <- FindVariableFeatures(fivemonth, selection.method = "vst", nfeatures = 2000)

#identify top 10 highly variable genes
top10 <- head(VariableFeatures(fivemonth), 10)

#plot variables features 
plot1variable <- VariableFeaturePlot(fivemonth)
plot2variable <- LabelPoints(plot = plot1variable, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1variable, plot2variable))

#Scaling the data
all.genes5mo <- rownames(fivemonth)
fivemonth <- ScaleData(fivemonth, features = all.genes5mo)

#performing linear pca
fivemonth <- RunPCA(fivemonth, features = VariableFeatures(object = fivemonth))
VizDimLoadings(fivemonth, dims = 1:2, reduction = "pca")
DimPlot(fivemonth, reduction = "pca")
DimHeatmap(fivemonth, dims = 1:15, cells = 500, balanced = TRUE)


#UMAP,  hange dimensions to test
#Clustering the cells
fivemonth <- FindNeighbors(fivemonth, dims = 1:12)
fivemonth <- FindClusters(fivemonth, resolution = 0.5)

#Running UMAP
fivemonth <- RunUMAP(fivemonth, dims = 1:12)
fivemonth <- RenameIdents(fivemonth, '0' = "Type IIb Myonuclei", '1' = "Type IIx Myonuclei", '2' = "Type IIb Myonuclei", '3' = "Type IIx Myonuclei", '4' = "FAPs", '5' = "Endothelial Cells", '6' = "Myotendinous Junction", '7' = "Smooth Muscle", '8' = "Satellite Cells", '9' = "Immune Cells", '10' = "Smooth Muscle", '11' = "Neuromuscular Junction", '12' = "Subcutaneous Fat")

DimPlot(fivemonth, reduction = "umap", label = TRUE)
fivemonth_muscle <- subset(fivemonth, idents = c("Type IIx Myonuclei", "Type IIb Myonuclei", "Satellite Cells", "Neuromuscular Junction", "Myotendinous Junction"))
fivemonth$CellType <- Idents(fivemonth)

#Finding Cluster markers
#fivemonth.markers <- FindAllMarkers(fivemonth, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#fivemonth_markers<- as.data.frame(fivemonth.markers) %>%
#  group_by(cluster) %>%
#  slice_max(n = 2, order_by = avg_log2FC)


#For 24 m/o mice
# 24-month data processing
# Get the data and create Seurat object
twentyfourmonth.data <- Read10X_h5("24-Month_filtered_feature_bc_matrix.h5", use.names = T)
twentyfourmonth <- CreateSeuratObject(counts = twentyfourmonth.data, project = "24-month", min.cells = 3, min.features = 200)
# Counting the percentage of reads that map to mitochondrial genome
twentyfourmonth [["percent.mt"]] <- PercentageFeatureSet(twentyfourmonth, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(twentyfourmonth, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Feature scatter plot
plot1 <- FeatureScatter(twentyfourmonth, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(twentyfourmonth, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# Filtering cells with nFeature RNA between 200 - 3200 and percentage of mitochondrial RNA < 5
twentyfourmonth <- subset(twentyfourmonth, subset = nFeature_RNA > 200 & nFeature_RNA < 3200 & percent.mt < 5)
# Normalize data with log transformation (default setting)
twentyfourmonth <- NormalizeData(twentyfourmonth, normalization.method = "LogNormalize", scale.factor = 10000)
# Identify highly variable features across cells (feature selection), return 2,000 features per dataset
twentyfourmonth <- FindVariableFeatures(twentyfourmonth, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(twentyfourmonth), 10)
# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(twentyfourmonth)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
# Scaling the data (linear transformation required for dimensional reduction techniques like PCA)
all.genes24mo <- rownames(twentyfourmonth)
twentyfourmonth <- ScaleData(twentyfourmonth, features = all.genes24mo)
# Perform linear dimensional reduction
twentyfourmonth <- RunPCA(twentyfourmonth, features = VariableFeatures(object = twentyfourmonth))
# Examine and visualize PCA results
print(twentyfourmonth[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(twentyfourmonth, dims = 1:2, reduction = "pca")
DimPlot(twentyfourmonth, reduction = "pca")
DimHeatmap(twentyfourmonth, dims = 1:15, cells = 500, balanced = TRUE)
# Cluster the cells
twentyfourmonth <- FindNeighbors(twentyfourmonth, dims = 1:15)
twentyfourmonth <- FindClusters(twentyfourmonth, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(twentyfourmonth), 5)
# Run non-linear dimensional reduction (UMAP/tSNE)
twentyfourmonth <- RunUMAP(twentyfourmonth, dims = 1:15)
DimPlot(twentyfourmonth, reduction = "umap",  label = TRUE)
twentyfourmonth <- RenameIdents(twentyfourmonth, "0" = "Type IIb Myonuclei", "1" = "Type IIx Myonuclei", "2" = "Type IIb Myonuclei #2", "3" = "Type IIx Myonuclei #2", "4" = "FAPs", "5" = "Smooth Muscle", "6" = "Endothelial Cells", "7" = "Myotendinous Junction", "8" = "Satellite Cells", "9" = "Immune Cells", "10" = "Smooth Muscle #2", "11" = "Endothelial Cells #2", "12" = "Neuromuscular Junction", "13" = "Schwann Cells", "14" = "Adipocytes")

DimPlot(twentyfourmonth, reduction = "umap", label = TRUE)
twentyfourmonth_muscle <- subset(twentyfourmonth, idents = c("Type IIx Myonuclei", "Type IIb Myonuclei", "Satellite Cells", "Neuromuscular Junction", "Myotendinous Junction", "Type IIb Myonuclei #2", "Type IIx Myonuclei #2"))

#integration od 5m/o and 24 m/o dataset
fivemonth_muscle@meta.data[["orig.ident"]] <- "A"
twentyfourmonth_muscle@meta.data[["orig.ident"]] <- "B"

integration.list.muscle <- list(fivemonth_muscle, twentyfourmonth_muscle)
for (i in 1:length(integration.list.muscle)) {
  integration.list.muscle[[i]] <- NormalizeData(integration.list.muscle[[i]], verbose = FALSE)
  integration.list.muscle[[i]] <- FindVariableFeatures(integration.list.muscle[[i]], selection.method = "vst", 
                                                       nfeatures = 2000, verbose = FALSE)
}

integration.anchors.muscle <- FindIntegrationAnchors(object.list = integration.list.muscle, dims = 1:30, verbose = FALSE)
integration.integrated.muscle <- IntegrateData(anchorset = integration.anchors.muscle, dims = 1:30)
integration.integrated.muscle <- ScaleData(integration.integrated.muscle, verbose = FALSE)
integration.integrated.muscle <- RunPCA(integration.integrated.muscle, npcs = 30, verbose = FALSE)
integration.integrated.muscle <- RunUMAP(integration.integrated.muscle, reduction = "pca", dims = 1:30)
DimPlot(integration.integrated.muscle, reduction = "umap", label = TRUE)
DimPlot(integration.integrated.muscle, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE) + NoLegend()

integration.integrated.muscle <- ScaleData(integration.integrated.muscle, verbose = FALSE)
muscle_markers <- FindAllMarkers(integration.integrated.muscle)
A <- subset(integration.integrated.muscle, orig.ident == "A")
B <- subset(integration.integrated.muscle, orig.ident == "B")

A_cells <- rownames(A)
B_cells <- rownames(B)

A_markers <- FindMarkers(integration.integrated.muscle, ident.1 = "A", ident.2 = NULL, group.by = "orig.ident", only.pos = TRUE)
B_markers <- FindMarkers(integration.integrated.muscle, ident.1 = "B", ident.2 = "A", group.by = "orig.ident", only.pos = TRUE)

top10 <- muscle_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.xlsx(A_markers, "5month positive markers.xlsx")
write.xlsx(B_markers, "24month positive markers.xlsx")

