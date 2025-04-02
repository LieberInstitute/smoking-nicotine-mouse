
# Script summary

## 01. SPEAQeasy 

Run [_SPEAQeasy_](https://doi.org/10.1186/s12859-021-04142-3) pipeline for quantification of expression features. 

* [01_make_manifest.R](01_SPEAQeasy/01_make_manifest.R): construct a manifest listing sample IDs and associated raw FASTQ files.

* [02_run_pipeline.sh](01_SPEAQeasy/02_run_pipeline.sh): invoke _SPEAQeasy_ to align data, quantify features, and produce [`RangedSummarizedExperiment`](https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html) objects for downstream analysis.



## 02. Build objects
This part of the code builds the necessary objects to analyze in downstream steps.

* [01_add_sequencing_lane.R](02_build_objects/01_add_sequencing_lane.R): extract flowcell information of the samples.

* [02_build_objects.R](02_build_objects/02_build_objects.R): explore, clean, correct, and format data. Here raw counts of genes, exons, and exon-exon junctions were log-normalized, and transcripts per million (TPM) of transcripts were log-scaled; lowly-expressed features were filtered out, and datasets separated by tissue and age. 



## 03. Exploratory Data Analysis (EDA)

* [01_StudyDesign.R](03_EDA/01_StudyDesign.R): compute the number of samples belonging to each pair/triad of sample variables.

* [02_QC.R](03_EDA/02_QC.R): evaluate and compare quality control (QC) metrics across sample groups; explore the relationships between QC metrics, and sample filtering by QC, followed by examination of QC metrics of samples that were kept and removed.

* [03_PCA_MDS.R](03_EDA/03_PCA_MDS.R): perform Principal Component Analysis (PCA) and Multidimensional Scaling Analysis (MDS) to explore sample-level variability in gene expression (at the 4 expression feature levels: genes, txs, exons, and jxns) and identify drivers of such variation. QC metrics of segregated samples in PCA plots were further analyzed and those samples were removed if turned out to be poor-quality. 

* [04_Expl_Var_partition.R](03_EDA/04_Expl_Var_partition.R): assess correlation between explanatory sample variables and explore their individual or joint contributions on the expression variation of each gene to select main contributing variables to be included in the models for differential expression. 



## 04. Differential Expression Analysis (DEA)
Separated in gene, transcript (tx), exon, and exon-exon junction (jx) level analyses. For each the following two scripts contain:

* At the [gene](04_DEA/Gene_analysis) level:
  * [01_Modeling.R](04_DEA/Gene_analysis/01_Modeling.R): perform differential expression analysis using [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html), separately for the 5 experimental groups: 
    * Nicotine exposure vs vehicle exposure in pup brain
    * Smoking exposure vs smoking control in pup brain
    * Nicotine administration vs vehicle administration in adult brain
    * Smoking exposure vs smoking control in adult brain
    * Smoking exposure vs smoking control in adult blood
    
    3 models were applied for each group: 
    * Naive model: gene expression is modeled by the following covariates: ~ `Group` + `plate` + `flowcell` + `QC metrics`
    * Adjusted model: gene expression is modeled by: ~ `Group` + [`Sex` (for pups) or `Pregnancy` (for adults)] + `plate` + `flowcell` + `QC metrics`
    * Interaction model: gene expression is modeled by: ~ [`Group*Sex` + `Group` + `Sex` (for pups) or `Group*Pregnancy` + `Group` + `Pregnancy` (for adults)] + `plate` + `flowcell` + `QC metrics`
    
    In addition, transcriptomic sex differences in pup brain in the nicotine and smoking experiments, were tested with differential expression analysis for `Sex` adjusting for `Group` and the same covariates, and testing differential expression for `Group` in female and male pups separately.
    
    `QC metrics` refer to any of `rRNA_rate`, `overallMapRate`, `totalAssignedGene`, `ERCCsumLogErr` and `mitoRate`. See [TableS18](../processed-data/SupplementaryTables/TableS18_sample_variable_dict.tsv) for the definition of these variables.

  * [02_Comparisons.R](04_DEA/Gene_analysis/02_Comparisons.R): 
  
    Create plots comparing the moderated *t*-stats of the genes in:                                                         
    **Analysis of blood vs brain biomarkers**
      * Smoking-exposed adult blood vs smoking/nicotine-exposed adult/pup brain 
      
          → search for brain genes (and pup brain transcripts) replicating in blood. 
          
    **Compare experiments**
      * Smoking-exposed pup brain vs nicotine-exposed pup brain
      * Smoking-exposed adult brain vs nicotine-exposed adult brain
      
    **Compare ages**     
      * Smoking-exposed pup brain vs smoking-exposed adult brain
      * Nicotine-exposed pup brain vs nicotine-exposed adult brain
      
    **Compare models for differential expression**       
      * Smoking-exposed pup brain under the naive model vs smoking-exposed pup brain under the fitted model
      * Nicotine-exposed pup brain under the naive model vs nicotine-exposed pup brain under the fitted model
      
    **Human vs mouse genes** 
      * Nicotine-exposed mouse adult/pup brain vs smoking-exposed human prenatal/adult brain 
      * Smoking-exposed mouse adult/pup brain vs smoking-exposed human prenatal/adult brain 
      * Smoking-exposed mouse adult blood vs smoking-exposed human prenatal/adult brain 
      
        → search for human brain genes that replicate in mouse brain/blood and vice versa.  
    
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

  * [01_Modeling.R](04_DEA/Tx_analysis/01_Modeling.R) for txs and [01_Modeling.R](04_DEA/Exon_analysis/01_Modeling.R) for exons: perform differential expression analysis applying only the fitted model for: 
    * Nicotine exposure vs vehicle exposure in pup brain
    * Smoking exposure vs smoking control in pup brain
    
    
  * [02_Comparisons.R](04_DEA/Tx_analysis/02_Comparisons.R) for txs and [02_Comparisons.R](04_DEA/Exon_analysis/02_Comparisons.R) for exons: 
  
    
    **Compare experiments**
    
     Plot the moderated *t*-stats of txs/exons in:
      * Smoking-exposed pup brain vs nicotine-exposed pup brain
      
    **Compare expression features**
      * Plot moderated *t*-stats of txs/exons in nicotine/smoking-exposed pup brain vs moderated *t*-stats of their genes in the same experiment
      
    For exons only:
      * Plot moderated *t*-stats of genes in nicotine/smoking-exposed pup brain vs mean of |*t*-stats of the gene - *t*-stats of the gene's exons|
      * Plot moderated *t*-stats of exons in nicotine/smoking-exposed pup brain vs |*t*-stats of the exon's gene - *t*-stats of the exon|
      
      
    Boxplots of expression lognorm counts of relevant genes and their txs/exons in the nicotine and smoking experiments. 
   
    Venn diagrams comparing:
      * Up/downregulated nicotine vs up/down smoking DE txs/exons
      * Up & downregulated nicotine vs up & down smoking DE txs/exons
      * Genes of up/downregulated nicotine/smoking DE txs/exons
      * DE txs'/exons' genes vs DEGs, with both features up/down in nicotine or smoking only, up/down in both experiments, up in nicotine and down in smoking, and up in smoking and down in nicotine
      * Compare DE txs' genes vs DE exons' genes vs DEGs in the same groups as above (up/down in either nic/smo only or in both experimets, or with different regulation directions in nic and smo)


    


* At the [exon-exon junction](04_DEA/Jx_analysis) level:

  * [01_Modeling.R](04_DEA/Jx_analysis/01_Modeling.R): differential expression analysis for junctions applying the fitted model for: 
    * Nicotine exposure vs vehicle exposure in pup brain
    * Smoking exposure vs smoking control in pup brain
    
  * [02_Comparisons.R](04_DEA/Jx_analysis/02_Comparisons.R): 
  
    Venn diagrams comparing:
      * DEGs vs DE txs' genes vs DE exons' genes vs all DE jxns' genes in nicotine and smoking
      * DEGs vs DE txs' genes vs DE exons' genes vs the genes associated with novel DE jxns in nicotine and smoking
      * DEGs vs DE txs' genes vs DE exons' genes vs the nearest, preceding, and following genes to novel DE jxns (without assigned gene), in nicotine and smoking
      * DEGs vs DE txs' genes vs DE exons' genes vs the genes of *AltStartEnd*, *InGen*, and *ExonSkip* DE jxns, in nicotine and smoking. See [caption of Supplementary Table 14](../processed-data/SupplementaryTables/) for an explanation of these exon-exon junction classes.
      

    Explore the mean expression and logFC of DEGs and non-DE genes with and without DE features in nicotine and smoking.   

 

## 05. Functional Enrichment Analysis (GO & KEGG terms)
Separated in gene, transcript (tx), and exon-level analyses.

* At the [gene](05_GO_KEGG/Gene_analysis) level:
  * [01_GO_KEGG_analyses.R](05_GO_KEGG/Gene_analysis/01_GO_KEGG_analyses.R): 
    * Overrepresentation analysis (ORA) of GO and KEGG terms among clusters of up and downregulated DEGs in nicotine and smoking pup brain (from fitted models only). 
    * Expression lognorm counts of the most significant genes within each cluster or in a certain enriched term (biological process (BP), molecular function (MF), cellular component (CC), or pathway (KEGG)) of interest were examined in boxplots. 
  
* At the [transcript](05_GO_KEGG/Tx_analysis) level:
  * [01_GO_KEGG_analyses.R](05_GO_KEGG/Tx_analysis/01_GO_KEGG_analyses.R): 
    * ORA of GO and KEGG terms among clusters of up and downregulated DE txs' genes in nicotine and smoking pup brain.
    * ORA of GO and KEGG terms among DE txs' genes not considered or non-DE at gene level in the nicotine and smoking pup brain.
    * ORA of GO and KEGG terms among DEGs and non-DE genes with DE exon and/or DE txs in the nicotine and smoking pup brain.
    * Create boxplots with lognorm counts of the most significant genes (with DE txs) in enriched terms of interest. 

* At the [exon](05_GO_KEGG/Exon_analysis) level:
  * [01_GO_KEGG_analyses.R](05_GO_KEGG/Exon_analysis/01_GO_KEGG_analyses.R): 
    * ORA of GO and KEGG terms among clusters of up and downregulated DE exons' genes in nicotine and smoking pup brain. 
    * Create boxplots with lognorm counts of the most significant genes (with DE exons) in enriched terms of interest. 



## 06. DGE visualization
Only at the gene level. 
 
* [01_Heatmap_DEG.R](06_Visualize_DEG/Gene_analysis/01_Heatmap_DEG.R): visualize gene expression patterns of all/up/downregulated pup brain DEGs in the nicotine and smoking experiments (from fitted models only) through heatmaps.


## 07. Novel DE junction gene annotation
* [01_Jxn_anno.R](07_Jxn_anno/01_Jxn_anno.R): here we explore jxn classes and their annotation in GENCODE M25; we quantify total and per-gene numbers of DE jxns from each class in nicotine and smoking. We search the nearest, preceding, and following genes of fully novel DE jxns (without associated gene) in nicotine and smoking pup brain. 


## 08. SRA submission
* [01_metadata.R](08_SRA/01_metadata.R) and [02_biosample.R](08_SRA/02_biosample.R): prepare metadata files as part of the upload of raw bulk RNA-seq FASTQ data to the Sequence Read. TODO 



