# mCIRCrna
The goal is to integrate circadian transcriptomes with muscle snRNA-seq data to identify age and cell type dependent circadian gene signatures that demonstrate tissue chronicity in muscle.


## Table of Contents

- [mCIRCrna](#team-repo-template)
    - [Background](#background)
    - [Hypothesis](#hypothesis)
    - [Aims](#aims)
    - [Data](#data)
    - [Method](#method)
        - Circadian Regulated Genes
        - Differential Expressions
        - Tools
    - [Results](#results)
        - Circadian Regulated Genes
        - Differential Expressions
    - [References](#references)
    - [Team Members](#team-members)

## Background

Skeletal muscle is an endocrine organ that composes 40-60% of the adult body mass with  metabolic and overall health span implications for the entire body. In aging it is well accepted that circadian rhythm declines and is dependent on transcription reprogramming. Because skeletal muscle is an endocrine organ consisting as a substantial portion of body mass, it poses as a promising organ for drug therapeutics, but its high cellular heterogeneity poses a challenge to safe druggable target. 

## Hypothesis

We hypothesize that we will be able to identify markers important to circadian rhythm and certain cell types if we cross- reference gene expression across age associated with circadian rhythm patterns from bulk transcriptomics to single cell muscle transcriptomics data.


## Aims

To identify circadian gene signatures that demonstrate tissue chronicity (2 fold expression changes over time) compared to age and cell type in muscle.

- Find genes that show chronicity (2 fold changes over time course ZT1-24) from CircAge RNA-seq transcriptomic database across age for young and old adult

- Cross reference signifcant genes that epress chronicity from CircAge RNA-seq transcriptomic to Myoatlas snRNA-seq muscle database for cell specific age related for young and old adult differential expression

## Data

- McCarthy JJ, Andrews JL, McDearmon EL, Campbell KS, Barber BK, Miller BH, Walker JR, Hogenesch JB, Takahashi JS, Esser KA. Identification of the circadian transcriptome in adult mouse skeletal muscle. Physiol Genomics. 2007 Sep 19;31(1):86-95. doi: 10.1152/physiolgenomics.00066.2007. Epub 2007 Jun 5. PMID: 17550994; PMCID: PMC6080860.

- CircAge: https://circaage.shinyapps.io/circaage/

- MyoAtlas: https://research.cchmc.org/myoatlas/

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

- Shufan Zhang | shufan0519@gmail.com | Team Leader  
- Lisa Shrestha | Member
- Van Nha Huynh | Member
- Kristen Coutinho | Member
- Russel Santos | Member
- Herminio Vazquez | Member

