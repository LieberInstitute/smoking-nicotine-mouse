smokingMouse 
================

## Citation
Code used to generate results for the writting: *"Modeling the effects of smoking and nicotine exposures on the developing brain"*
`TODO`

## Overview

This project consists of a differential expression analysis involving 4 data types: genes, exons, transcripts and junctions. The main goal of this study was to explore the effects of smoking and nicotine exposures on the developing brain of mice pups. As secondary objectives, this work evaluated the affected genes by each substance on adult brain in order to compare pup and adult results, and the effects of smoking on adult blood and brain to search for overlapping biomarkers in both tissues. 

## Study design
 
36 pregnant mice and 35 not pregnant female adults were administered nicotine (n=12), exposed to cigarette smoke (n=24) or controls (n=35) and RNA sequecing experiments were performed on frontal cortices of all the resultant 137 P0 pups and on frontal cortices (n=47) and blood (n=24) from the 71 adults, totaling 208 samples. Of the total pup samples, 19 were born to mice that were administered nicotine, 46 to mice exposed to smoking and the remainig 72 to control mice.

<img src="https://s3.us-west-2.amazonaws.com/secure.notion-static.com/618bb981-d4c7-4caa-8a1c-5af4cc0cfb83/Untitled.png?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Content-Sha256=UNSIGNED-PAYLOAD&X-Amz-Credential=AKIAT73L2G45EIPT3X45%2F20230131%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20230131T012754Z&X-Amz-Expires=86400&X-Amz-Signature=ee963a83bf64bc83fdd1ee01bf22428c0d0ef8d47e1dfb1a01e78f9bf31be85e&X-Amz-SignedHeaders=host&response-content-disposition=filename%3D%22Untitled.png%22&x-id=GetObject" width="600px" align="center" />

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
* *02_Comparisons.R*: Create plots comparing t-stats of different groups of features and Venn diagrams to quantify the number of common DE features between different experiments and/or regulation direction

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
  * `/dcl01/lieber/ajaffe/lab/smokingMouse_Indirects`
