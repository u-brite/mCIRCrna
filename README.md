# mCIRCrna
The goal is to integrate circadian transcriptomes with muscle snRNA-seq data to identify age and cell type dependent circadian gene signatures that demonstrate tissue chronicity in muscle.


## Table of Contents

- [mCIRCrna](#team-repo-template)
    - [Background](#Background)
    - [Data](#data)
    - [Method](#method)
        - Circadian Regulated Genes
        - Differential Expressions
        - Tools
    - [Results](#results) _Optional depending on project_
        - Circadian Regulated Genes
        - Differential Expressions
    - [References](#references)
    - [Team Members](#team-members)

## Background

Goal is to identify circadian gene signatures that demonstrate tissue chronicity (2 fold expression changes over time) compared to age and cell type in muscle.

Goal would be to use the CircAge RNA-seq transcriptomic database, and cross reference muscle transcriptomes with snRNA-seq muscle database. Genes that show chronicity (2 fold changes over time course ZT1-24) are compared to snRNA-seq profiles from Myoatlas. Expression profiles of dysregulated chronic genes (e.g. Per2) would then be compared in muscle cell populations (MTJ, FAPs, type II myonuclei, etc.). Conversely, any very dysregulated muscle genes across aging (P21, 24 months, 30 months TA) could then be cross referenced for chronicity using CircAge.

## Data

CircAge: https://circaage.shinyapps.io/circaage/

MyoAtlas: https://research.cchmc.org/myoatlas/

## Method

### Circadian regulated genes

### Differential Expression
Filtered feature bc files for each age group (5 month and 24 month) from snRNA Sequencing were downloaded from https://www.synapse.org/#!Synapse:syn21676145/files/ deposited by Millay DP et al. R (version 4.1.1) package Seurat (version 4.1.1) was used to explore the gene expressions across cell types and ages. For each age group, nuclei with less than 200 or greater than 3200 expressed features and the features that were expressed in less than three cells were excluded. Data were normalized logarithmically using the NormalizeData() function. The top 2000 features with the highest expression variability across the nuclei were identifided using FindVariableFeatures() function, which were then used in Principle Component Analysis (PCA) with RunPCA() function. For cell clustering and Uniform Manifold Approximation and Projection (UMAP) visualization, 12 and 15 components were used for 5 month and 24 month, respectively. The cell clustering and UMAP were done using FindNeighbors(), FindClusters(), and RunUMAP() functions. Cell types were assigned manually as did by Millay DP et al. The datasets were subset to select populations of interest, which are muscle cells. 
The 5 month and 24 month datasets were then integrated using FindIntegrationAnchors(), creating a new integrated Seurat object. PCA and UMAP visualization were done as described above. The function FindMarkers() was used to determine differential expression between populations of interest across ages and between cell populations. The FeaturePlot() function was used for visualization of the results. 


### Tools

## Results
:exclamation: _If your project yielded or intends to yield some novel analysis, please include them in your readme. It can be named something other than results as well._ :exclamation:

## References
    Ding, Haocheng, et al. Likelihood-based tests for detecting circadian rhythmicity and differential circadian patterns in transcriptomic applications. Briefings in Bioinformatics 22.6 (2021): bbab224.

## Team Members

Shufan Zhang | shufan0519@gmail.com | Team Leader  

