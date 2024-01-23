smokingMouse 
================

`Zenodo DOI:TODO`

Code used to generate results for the paper: *"Modeling the effects of smoking and nicotine exposures on the developing brain"* (Cite `TODO`)

## Citation
Original publication: *"Modeling the effects of smoking and nicotine exposures on the developing brain"* (Cite `DOI:TODO`)

## Overview

<style>
body {
text-align: justify}
</style>
This project consisted of a differential expression analysis involving 4 expression features: genes, exons, transcripts and exon-exon junctions. The main goal of this study was to explore the effects of prenatal smoking and nicotine exposures on the developing brain of mouse pups. As secondary objectives, this work evaluated the affected genes by each substance on adult brain in order to compare pup and adult results, and the effects of smoking exposure on adult blood and brain to search for overlapping biomarkers in both tissues. 

## Study design
 
<p align="center">
  <img src= "plots/03_EDA/01_StudyDesign/Study_design_fig.png" width="800" >
</p>

**Experimental design of the study.** **A)** 21 pregnant mice and 26 nonpregnant female adults were either administered nicotine (n=12), exposed to cigarette smoke (n=12), or used as controls (n=23; 11 nicotine controls and 12 smoking controls). A total of 137 pups were born to pregnant mice: 19 were born to mice that were administered nicotine, 46 to mice exposed to smoking, and the remaining 72 to control mice (23 to nicotine controls and 49 to smoking controls). Samples from frontal cortices of all P0 pups (n=137: 42 of nicotine and 95 of the smoking experiment) and adults (n=47: 23 of nicotine and 24 of the smoking experiment) were obtained, as well as blood samples from the smoking-exposed and smoking control adults (n=24), totaling 208 bulk RNA-seq samples. Number of donors and RNA-seq samples are indicated in the figure. **B)** RNA was extracted from such samples and RNA-seq experiments were performed, obtaining expression counts for genes, exons, transcripts and exon-exon junctions.


## Workflow

<p align="center">
  <img src= "plots/03_EDA/01_StudyDesign/Table_of_Analyses.png" width="1000" >
</p>

**Summary of analysis steps across gene expression feature levels**: 

* **1. Data processing**: counts of genes, exons, and exon-exon junctions were normalized to CPM and log2-transformed; transcript expression values were only log2-scaled since they were already in TPM. Lowly-expressed features were removed using the indicated functions and samples were separated by tissue and age in order to create subsets of the data for downstream analyses. 

* **2. Exploratory Data Analysis (EDA)**: QC metrics of the samples were examined and used to filter the poor quality ones. Sample level effects were explored through dimensionality reduction methods and segregated samples in PCA plots were removed from the datasets. Gene level effects were evaluated with analyses of variance partition. 

* **3. Differential Expression Analysis (DEA)**: with the relevant variables identified in the previous steps, the DEA was performed at the gene level for nicotine and smoking experiments in adult and pup brain samples, and for smoking in adult blood samples; DEA at the rest of the levels was performed for both exposures in pup brain only. DE signals of the genes in the different conditions, ages, tissues and species (human results from 1:[Semick et al. 2020](https://www.nature.com/articles/s41380-018-0223-1)) were contrasted, as well as the DE signals of exons and transcripts vs those of their genes. We also analyzed the mean expression of significant and non-significant genes with and without DE features. Then, all resultant DEGs and DE features (and their genes) were compared by direction of regulation (up or down) between and within experiments (nicotine/smoking); mouse DEGs were also compared against human genes associated with TUD from 2:[Toikumo et al. 2023](https://www.medrxiv.org/content/10.1101/2023.03.27.23287713v2). 

* **4. Functional Enrichment Analysis**: we obtained the GO & KEGG terms significantly enriched in our clusters of DEGs and genes of DE transcripts and exons.

* **5. DGE visualization**: the log2-normalized expression of DEGs was represented in heatmaps in order to distinguish the groups of up and downregulated genes.

* **6. Novel junction gene annotation**: for uncharacterized DE junctions with no annotated gene, their nearest, preceding and following genes were determined. See Materials and Methods for complete details. 

<p style="line-height:80%">
<font size="1.5"> 
Abbreviations: Jxn: junction; Tx(s): transcript(s); CPM: counts per million; TPM: transcripts per million; TMM: Trimmed Mean of M-Values; TMMwsp: TMM with singleton pairing; QC: quality control; PC: principal component; DEA: differential expression analysis; DE: differential expression/differentially expressed; FC: fold-change; FDR: false discovery rate; DEG: differentially expressed genes; TUD: tobacco use disorder; DGE: differential gene expression.
</font>
</p>

See [code](code/) for script summary. 


## Supplementary Tables

See [processed-data/SupplementaryTables](processed-data/SupplementaryTables/) to access and see a description of all supplementary tables generated in this study, including GitHub permalinks to the scripts in which they were created.


## Data access
`TODO`


## Internal 
* JHPCE locations:
  * `/dcs05/lieber/marmaypag/smokingMouseGonzalez2023_LIBD001`
  * `/dcl01/lieber/ajaffe/lab/smokingMouse_Indirects` (old location)
