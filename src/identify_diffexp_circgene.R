# set directory
getwd()
#setwd("/Users/lovely_shufan/Dropbox (Edison_Lab@UGA)/Hackin'Omics")
library(Seurat)
library(SeuratObject)

# input data
integration.integrated.muscle <- readRDS("intergration.integrated.muscle.rds")
df_clean2 <- read.csv("circadian_regulated_genes_by_circage.csv", header = TRUE) 

# differential expression between 2x and 2b muscle cells
# "0" = "Type IIb Myonuclei", "1" = "Type IIx Myonuclei", "2" = "Type IIb Myonuclei #2", 
# "3" = "Type IIx Myonuclei #2", "4" = "FAPs", "5" = "Smooth Muscle", "6" = "Endothelial Cells",
# "7" = "Myotendinous Junction", "8" = "Satellite Cells", "9" = "Immune Cells", "10" = "Smooth Muscle #2",
# "11" = "Endothelial Cells #2", "12" = "Neuromuscular Junction", "13" = "Schwann Cells", "14" = "Adipocytes"
library(ggplot2)
library(cowplot)
Idents(integration.integrated.muscle)
diffexp5month <- subset(integration.integrated.muscle, orig.ident == "A")
diffexp24month <- subset(integration.integrated.muscle, orig.ident == "B")

diffexp_5m_type2 <- FindMarkers(diffexp5month, ident.1 = "Type IIb Myonuclei", ident.2 = "Type IIx Myonuclei", verbose = FALSE)
diffexp_24m_type2 <- FindMarkers(diffexp24month, ident.1 = "Type IIb Myonuclei", ident.2 = "Type IIx Myonuclei", verbose = FALSE)
diffexp_24m_type2.2 <- FindMarkers(diffexp24month, ident.1 = "Type IIb Myonuclei #2", ident.2 = "Type IIx Myonuclei #2", verbose = FALSE)
diffexp_24m_type2b <- FindMarkers(diffexp24month, ident.1 = "Type IIb Myonuclei", ident.2 = "Type IIb Myonuclei #2", verbose = FALSE)
diffexp_24m_type2x <- FindMarkers(diffexp24month, ident.1 = "Type IIx Myonuclei", ident.2 = "Type IIx Myonuclei #2", verbose = FALSE)


# differential expression between Satellite vs myocyte
diffexp_5m_sattype2x <- FindMarkers(diffexp5month, ident.1 = "Satellite Cells", ident.2 = "Type IIx Myonuclei", verbose = FALSE)
diffexp_5m_sattype2b <- FindMarkers(diffexp5month, ident.1 = "Satellite Cells", ident.2 = "Type IIb Myonuclei", verbose = FALSE)
diffexp_24m_sattype2b <- FindMarkers(diffexp24month, ident.1 = "Satellite Cells", ident.2 = "Type IIb Myonuclei", verbose = FALSE)
diffexp_24m_sattype2x <- FindMarkers(diffexp24month, ident.1 = "Satellite Cells", ident.2 = "Type IIx Myonuclei", verbose = FALSE)
diffexp_24m_sattype2b.2 <- FindMarkers(diffexp24month, ident.1 = "Satellite Cells", ident.2 = "Type IIb Myonuclei #2", verbose = FALSE)
diffexp_24m_sattype2x.2 <- FindMarkers(diffexp24month, ident.1 = "Satellite Cells", ident.2 = "Type IIx Myonuclei #2", verbose = FALSE)

gene.names = c()
gene.names = c(gene.names, rownames(diffexp_5m_type2), rownames(diffexp_24m_type2), rownames(diffexp_24m_type2.2))
gene.names = c(gene.names, rownames(diffexp_24m_type2b), rownames(diffexp_24m_type2x), rownames(diffexp_5m_sattype2x))
gene.names = c(gene.names, rownames(diffexp_5m_sattype2b), rownames(diffexp_24m_sattype2b), rownames(diffexp_24m_sattype2x))

gene.names = c(gene.names, rownames(diffexp_24m_sattype2b.2), rownames(diffexp_24m_sattype2x.2))
unique.gene.names = unique(gene.names)
diffexp.circgene.celltype = intersect(unique.gene.names, unique(df_clean2$Gene)) 

# differential expression by timepoints
integration.integrated.muscle$celltype.age <- paste(Idents(integration.integrated.muscle), integration.integrated.muscle$orig.ident, sep = "_")
integration.integrated.muscle$celltype <- Idents(integration.integrated.muscle)
Idents(integration.integrated.muscle) <- "celltype.age"
typeIIb.age <- FindMarkers(integration.integrated.muscle, ident.1 = "Type IIb Myonuclei_A", ident.2 = "Type IIb Myonuclei_B", verbose = FALSE)
typeIIx.age <- FindMarkers(integration.integrated.muscle, ident.1 = "Type IIx Myonuclei_A", ident.2 = "Type IIx Myonuclei_B", verbose = FALSE)
satellitte.age <- FindMarkers(integration.integrated.muscle, ident.1 = "Satellite Cells_A", ident.2 = "Satellite Cells_B", verbose = FALSE)
myotendinous.age <- FindMarkers(integration.integrated.muscle, ident.1 = "Myotendinous Junction_A", ident.2 = "Myotendinous Junction_B", verbose = FALSE)
neuromuscular.age <- FindMarkers(integration.integrated.muscle, ident.1 = "Neuromuscular Junction_A", ident.2 = "Neuromuscular Junction_B", verbose = FALSE)


gene.names.age = c(rownames(typeIIb.age), rownames(typeIIx.age), rownames(myotendinous.age), rownames(neuromuscular.age))
diffexp.circgene.age = intersect(unique(gene.names.age), unique(df_clean2$Gene)) 
diffexp.circgene.celltype.age = intersect(diffexp.circgene.age, diffexp.circgene.celltype)

