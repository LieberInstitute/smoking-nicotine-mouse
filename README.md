smokingMouse 
================

## Citation
Code used to generate results for the paper: *"Modeling the effects of smoking and nicotine exposures on the developing brain"*
`DOI:TODO`

## Overview

This project consisted of a differential expression analysis involving 4 expression features: genes, exons, transcripts and exon-exon junctions. The main goal of this study was to explore the effects of prenatal smoking and nicotine exposures on the developing brain of mouse pups. As secondary objectives, this work evaluated the affected genes by each substance on adult brain in order to compare pup and adult results, and the effects of smoking exposure on adult blood and brain to search for overlapping biomarkers in both tissues. 

## Study design
 
<p align="center">
  <img src= "plots/03_EDA/01_StudyDesign/Study_design_fig.png" width="800" >
</p>

**Experimental design of the study.** **A)** 21 pregnant mice and 26 nonpregnant female adults were either administered nicotine (n=12), exposed to cigarette smoke (n=12), or used as controls (n=23; 11 nicotine controls and 12 smoking controls). A total of 137 pups were born to pregnant mice: 19 were born to mice that were administered nicotine, 46 to mice exposed to smoking, and the remaining 72 to control mice (23 to nicotine controls and 49 to smoking controls). Samples from frontal cortices of all P0 pups (n=137: 42 of nicotine and 95 of the smoking experiment) and adults (n=47: 23 of nicotine and 24 of the smoking experiment) were obtained, as well as blood samples from the smoking-exposed and smoking control adults (n=24), totaling 208 bulk RNA-seq samples. Number of donors and RNA-seq samples are indicated in the figure. **B)** RNA was extracted from such samples and RNA-seq experiments were performed, obtaining expression counts for genes, exons, transcripts and exon-exon junctions.


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


## Data access
`TODO`


## Internal 
`TODO`
* JHPCE locations:
  * `/dcs05/lieber/marmaypag/smokingMouseGonzalez2023_LIBD001`
  * (old location) `/dcl01/lieber/ajaffe/lab/smokingMouse_Indirects`
