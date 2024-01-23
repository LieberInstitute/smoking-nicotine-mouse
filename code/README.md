
## Script summary

### 1. Build objects
This initial part of the code builds the necessary objects to analyze in the EDA. 
* *01_add_sequencing_lane.R*: Add flowcell information of the samples
* *02_build_objects.R*: Explore datasets, data cleaning, normalization and separation, and feature filtering

### 2. Exploratory Data Analysis
* *01_StudyDesign.R*: Compute the number of samples belonging to each pair/trio of sample phenotypes
* *02_QC.R*: Explore the relatioships between QC metrics and sample phenotypes, plot QC metrics per sample and sample filtering by QC
* *03_PCA_MDS.R*: Explore sample effects: Dimensionality reduction
* *04_Expl_Var_partition.R*: Explore gene level effects: Explanatory variables and Variance partition

### 3. Differential Expression Analysis
Separated in gene, tx, exon and jx level analyses
* *01_Modeling.R*: Perform DEA 
* *02_Comparisons.R*: Create plots comparing t-stats of different groups of features and Venn diagrams to quantify the number of common DE features between different experiments and/or regulation directions. Analysis of blood vs brain biomarkers and comparisons of human vs mouse genes are within the gene analysis. 

### 4. GO & KEGG analyses
Separated in gene, tx and exon level analyses
* *01_GO_KEGG_analyses.R*: Gene Ontology and KEGG enrichment analyses at the gene level

### 5. DEG visualization
Only at the gene level
* *01_Heatmap_DEG.R*: Create heatmaps of DEG 

### 6. Junction annotation
* *01_Jxn_anno.R*: Obtain novel DE jxns and explore potential novel isoforms
