
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
Separated in gene, transcript (Tx), exon and exon-exon junction (Jx) level analyses. For each the following two scripts contain:

* At the [gene](04_DEA/Gene_analysis) level:
  * [01_Modeling.R](04_DEA/Gene_analysis/01_Modeling.R): Perform differential expression analysis using [`limma`](https://bioconductor.org/packages/release/bioc/html/limma.html), separately for the 5 experimental groups: 
    * Smoking-exposed vs smoking control adult blood
    * Nicotine-exposed vs nicotine control adult brain
    * Smoking-exposed vs smoking control adult brain
    * Nicotine-exposed vs nicotine control pup brain
    * Smoking-exposed vs smoking control pup brain  
    
    3 models were applied for each group: 
    * Naive model: gene expression is modeled by the following covariates: ~ `Group` + `plate` + `flowcell` + `QC metrics`
    * Fitted model: gene expression is modeled by: ~ `Group` + [`Sex` (for pups) or `Pregnancy` (for adults)] + `plate` + `flowcell` + `QC metrics`
    * Interaction model: gene expression is modeled by: ~ [`Group*Sex` (for pups) or `Group*Pregancy` (for adults)] + `plate` + `flowcell` + `QC metrics`
    

  * [02_Comparisons.R](04_DEA/Gene_analysis/02_Comparisons.R): 
  
    Create plots comparing the moderated *t*-stats of the genes in:
   
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
      * Nicotine-exposed mouse adult/pup brain vs smoking-exposed human prenatal/adult brain 
      * Smoking-exposed mouse adult/pup brain vs smoking-exposed human prenatal/adult brain 
      * Smoking-exposed mouse adult blood vs smoking-exposed human prenatal/adult brain 
      
        → search for human brain genes that replicate in mouse brain/blood and vice versa

     Venn diagrams to quantify the number of:
      * Nicotine/smoking DEGs in the naive vs fitted model 
      * Naive/fitted model DEGs in nicotine vs smoking
      * Up/downregulated nicotine/smoking DEGs in the naive vs fitted model 
      * Up/downregulated nicotine DEGs vs up/down smoking DEGs (from either naive or fitted model or from only one model)
      * Nicotine naive & fitted DEGs vs smoking naive & fitted DEGs
      * Up & downregulated nicotine DEGs vs up & down smoking DEGs (from the fitted model only)
      * All nicotine vs smoking DEGs (from any model)
      * Mouse genes replicating in human
      * Tobacco Use Disorder (TUD)-associated human genes vs mouse pup DEGs 



* At the [transcript](04_DEA/Tx_analysis) and [exon](04_DEA/Exon_analysis) levels:

  * [01_Modeling.R](04_DEA/Tx_analysis/01_Modeling.R) for txs and [01_Modeling.R](04_DEA/Exon_analysis/01_Modeling.R) for exons: Perform differential expression analysis applying only the fitted model for: 
    * Nicotine-exposed vs nicotine control pup brain
    * Smoking-exposed vs smoking control pup brain  
    
    
  * [02_Comparisons.R](04_DEA/Tx_analysis/02_Comparisons.R) for txs and [02_Comparisons.R](04_DEA/Exon_analysis/02_Comparisons.R) for exons: 
  
    Plot the moderated *t*-stats of txs/exons in:
    
    **Compare experiments**
      * Smoking-exposed pup brain vs nicotine-exposed pup brain
      
    **Compare expression features**
      * Moderated *t*-stats of txs/exons in nicotine/smoking-exposed pup brain vs moderated *t*-stats of their genes in the same experiment
      
    For exons only:
      * Moderated *t*-stats of genes in nicotine/smoking-exposed pup brain vs mean of |*t*-stats of the gene - *t*-stats of the gene's exons|
      * Moderated *t*-stats of exons in nicotine/smoking-exposed pup brain vs |*t*-stats of the exon's gene - *t*-stats of the exon|
      
      
    Boxplots of expression lognorm counts of relevant genes and their txs/exons in the nicotine and smoking experiments. 
   
    Venn diagrams comparing:
      * Up/downregulated nicotine vs up/down smoking DE txs/exons
      * Up & downregulated nicotine vs up & down smoking DE txs/exons
      * Genes of up/downregulated nicotine/smoking DE txs/exons
      * DE txs'/exons' genes vs DEGs, with both features up/down in nicotine or smoking only, up/down in both experiments, up in nicotine and down in smoking, and up in smoking and down in nicotine
      * Compare DE txs' genes vs DE exons' genes vs DEGs in the same groups as above (up/down in either nic/smo only or in both experimets, or with different regulation directions in nic and smo)


    


* At the [exon-exon junction](04_DEA/Jx_analysis) level:

  * [01_Modeling.R](04_DEA/Jx_analysis/01_Modeling.R): Differential expression analysis for junctions applying the fitted model for: 
    * Nicotine-exposed vs nicotine control pup brain
    * Smoking-exposed vs smoking control pup brain 
    
  * [02_Comparisons.R](04_DEA/Jx_analysis/02_Comparisons.R): 
  
    Venn diagrams comparing:
      * Up/downregulated nicotine vs up/down smoking DE txs/exons    
    

 

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
