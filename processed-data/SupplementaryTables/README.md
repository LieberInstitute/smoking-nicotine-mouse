# Supplementary Tables

## [Supplementary Table 1](TableS1_DEGs_brain_pup_nicotine.tsv)

**DEGs in the nicotine-exposed pup brain.**

Metadata of DEGs in the nicotine pup brain and their logFC, moderated *t*-stats, *p*-value and adjusted *p*-value for DE. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.

## [Supplementary Table 2](TableS2_DEGs_brain_pup_smoking.tsv)
**DEGs in the smoking-exposed pup brain.** 

Same as in [**Table S1**](TableS1_DEGs_brain_pup_nicotine.tsv) but for DEGs in the smoking pup brain.

## [Supplementary Table 3](TableS3_DEA_results_all_genes.tsv)
**Differential gene expression results for the complete gene dataset.** 

Gene-level metadata and the logFC, moderated *t*-stats, *p*-value and adjusted *p*-value for DE of each gene in the 5 experimental groups: smoking-exposed adult blood, and smoking/nicotine-exposed adult/pup brain. Also included are the replication results of the genes in mouse blood. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.

## [Supplementary Table 4](TableS4_DE_txs_brain_pup_nicotine.tsv)
**Differentially expressed transcripts in the nicotine-exposed pup brain.**

Metadata of DE transcripts in the nicotine pup brain and their logFC, moderated *t*-stats, *p*-value and adjusted *p*-value. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.

## [Supplementary Table 5](TableS5_DE_txs_brain_pup_smoking.tsv)
**Differentially expressed transcripts in the smoking-exposed pup brain.**

Same as in [**Table S4**](TableS4_DE_txs_brain_pup_nicotine.tsv) but for DE transcripts in the smoking pup brain.

## [Supplementary Table 6](TableS6_DE_txs_vs_genes_nic.tsv)
**Differential expression of transcripts vs genes for the nicotine experiment in pup brain.**

DE statistics (logFC, moderated *t*-stats, *p*-value and adjusted *p*-value) for transcripts and their respective genes in nicotine-exposed pup brain samples, and if transcripts and genes were both or solely DE. Only transcripts of genes present in the gene dataset are shown. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.

## [Supplementary Table 7](TableS7_DE_txs_vs_genes_smo.tsv)
**Differential expression of transcripts vs genes for the smoking experiment in pup brain.**

Same as in [**TableS6**](TableS6_DE_txs_vs_genes_nic.tsv) but for the smoking experiment.

## [Supplementary Table 8](TableS8_DE_exons_brain_pup_nicotine.tsv)
**Differentially expressed exons in the nicotine-exposed pup brain.**

Metadata of DE exons in the nicotine pup brain and their logFC, moderated *t*-stats, *p*-value and adjusted *p*-value. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.

## [Supplementary Table 9](TableS9_DE_exons_brain_pup_smoking.tsv) 
**Differentially expressed exons in the smoking-exposed pup brain.**

Same as in [**Table S8**](TableS8_DE_exons_brain_pup_nicotine.tsv) but for DE exons in the smoking pup brain.

## [Supplementary Table 10](TableS10_DE_exons_vs_genes_nic.tsv)
**Differential expression of exons vs genes for the nicotine experiment in pup brain.**

logFC, moderated *t*-stats, *p*-value and adjusted *p*-value for exons and their respective genes in nicotine-exposed pup brain samples, as well as if exons and genes were both DE or not. Only exons of genes present in the gene dataset are shown. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.

## [Supplementary Table 11](TableS11_DE_exons_vs_genes_smo.tsv)
**Differential expression of exons vs genes for the smoking experiment in pup brain.**

Same as in [**Table S10**](TableS10_DE_exons_vs_genes_nic.tsv) but for the smoking experiment.

## [Supplementary Table 12](TableS12_DE_jxns_brain_pup_nicotine.tsv) 
**Differentially expressed exon-exon junctions in the nicotine-exposed pup brain.**

Metadata of DE exon-exon junctions in the nicotine pup brain, including for each: 

  * if both its donor and acceptor sites together are known and annotated in GENCODE M25 (`inGencode` variable); 
  * if its donor and acceptor sites are individually annotated in GENCODE M25 (`inGencodeStart` and `inGencodeEnd` variables, respectively); 
  * the junction `Class`: *Novel* (if both start and end sites are unknown, also known as fully novel junctions), *InGen* (already annotated in GENCODE M25), *AltStartEnd* (if it has only one known site), or *ExonSkip* (with sites from non-successive exons, both known individually but not together), and
  * if it is a fusion junction, meaning that it connects exons from different genes (`isFusion` variable). 

  Their logFC, moderated *t*-stats, *p*-value and adjusted *p*-value are provided. These statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.


## [Supplementary Table 13](TableS13_DE_jxns_brain_pup_smoking.tsv)
**Differentially expressed exon-exon junctions in the smoking-exposed pup brain.**

Same as in [**Table S12**](TableS12_DE_jxns_brain_pup_nicotine.tsv) but for DE exon-exon junctions in the smoking pup brain.


## [Supplementary Table 14](TableS14_novel_DE_jxns_genes.tsv)
**Associated genes of fully novel DE exon-exon junctions in pup brain.**

Nearest, following and preceding genes of the fully novel DE exon-exon junctions without assigned gene in the nicotine and smoking-exposed pup brains.


## [Supplementary Table 15](TableS15_human_vs_mice_results.tsv)
**Differential gene expression results for gene pairs of mouse-human homologs.**

The logFC, moderated *t*-stats, *p*-value and adjusted *p*-value for the human gene in smoking-exposed prenatal and adult human brain, and for the mouse gene in the 5 experimental mice groups (as in [**Table S3**](TableS3_DEA_results_all_genes.tsv)), are presented. Only mouse genes with human homolog(s) present in the human dataset from [Semick et al. 2020](https://www.nature.com/articles/s41380-018-0223-1) are considered. The DE statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.


## [Supplementary Table 16](TableS16_TUD_human_genes_vs_mouseDEGs.tsv)
**Mouse DEGs in pup brain with human homologs TUD-associated.**

Pup brain DEGs in the nicotine and smoking conditions with human homologous genes that were the nearest genes of genome-wide significant (GWS) lead SNPs in loci associated with Tobacco Use Disorder (TUD), obtained from a multi-ancestral GWAS meta-analysis of TUD cases and controls in individuals from 8 cohorts (including UKBB), with European (EUR), African American (AA) and Latin American (LA) ancestry (*TUD-multi+UKBB* dataset), and from a within-ancestry GWAS meta-analysis in EUR individuals from 5 cohorts, including UKBB data (*TUD-EUR+UKBB* dataset). As well as human genes significantly associated with TUD in EUR individuals (*TUD-EUR-MAGMA* dataset); neurobiologically relevant target human genes associated with TUD (*TUD-EUR-H-MAGMA* dataset), especially expressed in prenatal (*TUD-EUR-H-MAGMA-prenatal* dataset) and adult brain (*TUD-EUR-H-MAGMA-adult* dataset); human genes whose expression is predicted to be affected by the EUR-SNPs across multiple brain regions (*TUD-EUR-S-MultiXcan* dataset) and in specific brain regions (*TUD-EUR-S-PrediXcan* dataset), including the frontal cortex (*TUD-EUR-S-PrediXcan-FC* dataset). See more details of these TUD-associated human genes in the original publication ([Toikumo et al. 2023](https://www.medrxiv.org/content/10.1101/2023.03.27.23287713v2)). 

## [Supplementary Table 17](TableS17_sample_variable_dict.tsv)
**Dictionary of sample variables.**

Description of the sample variables used throughout this project. 
 
 
## [Supplementary Table 18](TableS18_study_design.tsv)
**Study design.**

Number of samples from each pair of sample-level variables. See Table S17 for sample variable description


## [Supplementary Table 19](TableS19_study_design_afterSampleFilters.tsv)
**Samples used for downstream analyses.** 

Number of samples from each pair of sample-level variables after sample filtering based on QC metrics and PCA plots. See [**Table S17**](TableS17_sample_variable_dict.tsv) for sample variable description.



