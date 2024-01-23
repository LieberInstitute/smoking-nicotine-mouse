
# Script summary

## 01. SPEAQeasy 

Run [_SPEAQeasy_](https://doi.org/10.1186/s12859-021-04142-3) pipeline for quantification of expression features. 

* [01_make_manifest.R](01_SPEAQeasy/01_make_manifest.R): construct a manifest listing sample IDs and associated raw FASTQ files.

* [02_run_pipeline.sh](01_SPEAQeasy/02_run_pipeline.sh): invoke _SPEAQeasy_ to align data, quantify features, and produce [`RangedSummarizedExperiment`](https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html) objects for downstream analysis.



## 02. Build objects
This part of the code builds the necessary objects to analyze in downstream steps.

* [01_add_sequencing_lane.R](02_build_objects/01_add_sequencing_lane.R): Extract flowcell information of the samples.

* [02_build_objects.R](02_build_objects/02_build_objects.R): Explore, clean, correct and format data. Here raw counts are log-normalized, lowly-expressed features filtered out and datasets separated by tissue and age. 



## 03. Exploratory Data Analysis (EDA)

* [01_StudyDesign.R](03_EDA/01_StudyDesign.R): Compute the number of samples belonging to each pair/triad of sample phenotypes.

* [02_QC.R](03_EDA/02_QC.R): Evaluate and compare quality control (QC) metrics of sample groups; explore the relationships between QC metrics, and sample filtering by QC, followed by examination of QC metrics of samples that were kept and removed.

* [03_PCA_MDS.R](03_EDA/03_PCA_MDS.R): Perform Principal Component Analysis (PCA) and Multidimensional Scaling Analysis (MDS) to explore sample-level variability in gene expression and identify drivers of such variation. QC metrics of segregated samples in PCA plots were further analyzed and those samples were removed if turned out to be poor-quality. 

* [04_Expl_Var_partition.R](03_EDA/04_Expl_Var_partition.R): Assess correlation between explanatory sample variables and explore their individual or joint contributions on the expression variation of each gene to select main contributing variables to be included in the models for differential expression. 



## 04. Differential Expression Analysis (DEA)
Separated in gene, transcript (Tx), [exon](04_DEA/Exon_analysis) and [exon-exon junction](04_DEA/Jx_analysis) (Jx) level analyses. For each the following two scripts contain:

* At the [gene](04_DEA/Gene_analysis) level:
  * [01_Modeling.R](04_DEA/Gene_analysis/01_Modeling.R): Perform differential expression analysis using [`limma`](https://bioconductor.org/packages/release/bioc/html/limma.html) separately for the 5 experimental groups: 
    * Smoking-exposed vs smoking controls in adult blood
    * Nicotine-exposed vs nicotine controls in adult brain
    * Smoking-exposed vs smoking controls in adult brain
    * Nicotine-exposed vs nicotine controls in pup brain
    * Smoking-exposed vs smoking controls in pup brain  
    
    3 models were applied for each group: 
    * Naive model: gene expression is modeled by the following covariates: ~ `Group` + `plate` + `flowcell` + `QC metrics`
    * Fitted model: gene expression is modeled by: ~ `Group` + [`Sex` (for pups) or `Pregnancy` (for adults)] + `plate` + `flowcell` + `QC metrics`
    * Interaction model: gene expression is modeled by: ~ [`Group*Sex` (for pups) or `Group*Pregancy` (for adults)] + `plate` + `flowcell` + `QC metrics`
    

  * [02_Comparisons.R](04_DEA/Gene_analysis/02_Comparisons.R): Create plots comparing the moderated *t*-stats of genes in:
   
    **Analysis of blood vs brain biomarkers**
      * Smoking-exposed adult blood vs smoking/nicotine-exposed adult/pup brain 
      
          → search for brain genes (and pup brain transcripts) replicating in blood
      
    **Compare experiments**
      * Smoking-exposed pup brain vs nicotine-exposed pup brain
      * Smoking-exposed adult brain vs nicotine-exposed adult brain
      
    **Compare ages**     
      * Smoking-exposed pup brain vs smoking-exposed adult brain
      * Nicotine-exposed pup brain vs nicotine-exposed adult brain
      
    **Compare models for differential expression**       
      * Smoking-exposed pup brain in naive model vs smoking-exposed pup brain in fitted model
      * Nicotine-exposed pup brain in naive model vs nicotine-exposed pup brain in fitted model
      
    **Human vs mouse genes** 
      * Nicotine-exposed adult/pup brain vs smoking-exposed human prenatal/adult brain 
      * Smoking-exposed adult/pup brain vs smoking-exposed human prenatal/adult brain 
      * Smoking-exposed adult blood vs smoking-exposed human prenatal/adult brain 
      
        → search for human brain genes that replicate in mouse brain/blood and viceversa

     Venn diagrams to quantify the number of:
      * Nicotine/smoking DEGs in the naive vs fitted model 
      * Naive/fitted model DEGs in nicotine vs smoking
      * Up/downregulated nicotine/smoking DEGs in the naive vs fitted model 
      * Up/downregulated nicotine DEGs vs up/down smoking DEGs (from either naive or fitted model or from only one model)
      * Nicotine/smoking naive/fitted DEGs 
      * Up/downregulated nicotine/smoking DEGs (from the fitted model only)
      * All nicotine vs smoking DEGs (from any model)
      * Mouse replicating genes in human
      * Tobacco Use Disorder (TUD) associated human genes vs mouse pup DEGs 

    

* At the [transcript](04_DEA/Tx_analysis) level:



## 05. Functionl Enrichment Analysis (GO & KEGG)
Separated in gene, transcript (Tx), exon and exon-exon junction (Jx) level analyses.

* *01_GO_KEGG_analyses.R*: Gene Ontology and KEGG enrichment analyses at the gene level



## 06. DGE visualization
Only at the gene level. 

* *01_Heatmap_DEG.R*: Create heatmaps of DEG 



## 07. Novel DE junction gene annotation
* *01_Jxn_anno.R*: Obtain novel DE jxns and explore potential novel isoforms



# Code organization

Files are organized following the structure from [LieberInstitute/template_project](https://github.com/LieberInstitute/template_project). Scripts include the  R session information with details about version numbers of the packages we used.
