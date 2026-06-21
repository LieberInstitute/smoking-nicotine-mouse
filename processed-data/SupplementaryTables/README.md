# Supplementary Tables

## [Supplementary Table 1](TableS1_study_design.tsv)
**Study design.**

Number of samples from each pair of sample-level variables. See [**Table S19**](TableS19_sample_variable_dict.tsv) for sample variable description.



## [Supplementary Table 2](TableS2_study_design_afterSampleFilters.tsv)
**Samples used for downstream analyses.** 

Number of samples from each pair of sample-level variables after sample filtering based on QC metrics and PCA plots. See [**Table S19**](TableS19_sample_variable_dict.tsv) for sample variable description.



## [Supplementary Table 3](TableS3_DEGs_brain_pup_nicotine.tsv)

**DEGs in the nicotine-exposed pup brain.**

Metadata of DEGs in the nicotine pup brain and their logFC, moderated *t*-stats, *p*-value, and adjusted *p*-value for DE. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names. Data from this table were used to create the volcano plot in **Figure 2A**.

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Gene_analysis/01_Modeling.R#L468).



## [Supplementary Table 4](TableS4_DEGs_brain_pup_smoking.tsv)
**DEGs in the smoking-exposed pup brain.** 

Same as in [**Table S3**](TableS3_DEGs_brain_pup_nicotine.tsv) but for DEGs in the smoking pup brain. Data from this table were used to create the volcano plot in **Figure 2B**.

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Gene_analysis/01_Modeling.R#L524).



## [Supplementary Table 5](TableS5_DEA_results_all_genes.tsv)
**Differential gene expression results for the complete gene dataset.** 

Gene-level metadata and the logFC, moderated *t*-stats, *p*-value, and adjusted *p*-value for DE of each gene in the 5 experimental groups: nicotine vs vehicle exposure in pup brain, smoking exposure vs control in pup brain, nicotine vs vehicle administration in adult brain, smoking exposure vs control in adult brain, and smoking exposure vs control in adult blood. Also included are the replication results of the genes in mouse blood. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names. Scatter plots of *t*-stats in **Figure 2C**, **Figure 4** and **Figure S18A-B** can be reproduced with the data provided in this table.

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Gene_analysis/01_Modeling.R#L561).



## [Supplementary Table 6](TableS6_DE_txs_brain_pup_nicotine.tsv)
**Differentially expressed transcripts in the nicotine-exposed pup brain.**

Metadata of DE transcripts in the nicotine pup brain and their logFC, moderated *t*-stats, *p*-value, and adjusted *p*-value. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Tx_analysis/01_Modeling.R#L290).



## [Supplementary Table 7](TableS7_DE_txs_brain_pup_smoking.tsv)
**Differentially expressed transcripts in the smoking-exposed pup brain.**

Same as in [**Table S6**](TableS6_DE_txs_brain_pup_nicotine.tsv) but for DE transcripts in the smoking pup brain.

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Tx_analysis/01_Modeling.R#L313).



## [Supplementary Table 8](TableS8_DE_txs_vs_genes_nic.tsv)
**Differential expression of transcripts vs genes for the nicotine experiment in pup brain.**

DE statistics (logFC, moderated *t*-stats, *p*-value, and adjusted *p*-value) for transcripts and their respective genes for nicotine exposure in pup brain, and if transcripts and genes were both or solely DE. Only transcripts of genes present in the gene dataset are shown. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names. Scatter plot of nicotine *t*-stats in **Figure 3A** can be reproduced with the data contained in this table.  

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Tx_analysis/02_Comparisons.R#L514).



## [Supplementary Table 9](TableS9_DE_txs_vs_genes_smo.tsv)
**Differential expression of transcripts vs genes for the smoking experiment in pup brain.**

Same as in [**TableS8**](TableS8_DE_txs_vs_genes_nic.tsv) but for the smoking experiment. Scatter plot of smoking *t*-stats in **Figure 3A** can be reproduced with the data contained in this table. 

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Tx_analysis/02_Comparisons.R#L524).



## [Supplementary Table 10](TableS10_DE_exons_brain_pup_nicotine.tsv)
**Differentially expressed exons in the nicotine-exposed pup brain.**

Metadata of DE exons in the nicotine pup brain and their logFC, moderated *t*-stats, *p*-value, and adjusted *p*-value. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Exon_analysis/01_Modeling.R#L304).



## [Supplementary Table 11](TableS11_DE_exons_brain_pup_smoking.tsv) 
**Differentially expressed exons in the smoking-exposed pup brain.**

Same as in [**Table S10**](TableS10_DE_exons_brain_pup_nicotine.tsv) but for DE exons in the smoking pup brain.

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Exon_analysis/01_Modeling.R#L327).



## [Supplementary Table 12](TableS12_DE_exons_vs_genes_nic.tsv)
**Differential expression of exons vs genes for the nicotine experiment in pup brain.**

DE statistics (logFC, moderated *t*-stats, *p*-value, and adjusted *p*-value) for exons and their respective genes for nicotine exposure in pup brain, as well as if exons and genes were both DE or not. Only exons of genes present in the gene dataset are shown. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names. Scatter plot of nicotine *t*-stats in **Figure 3B** can be reproduced with the data contained in this table. 

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Exon_analysis/02_Comparisons.R#L296).



## [Supplementary Table 13](TableS13_DE_exons_vs_genes_smo.tsv)
**Differential expression of exons vs genes for the smoking experiment in pup brain.**

Same as in [**Table S12**](TableS12_DE_exons_vs_genes_nic.tsv) but for the smoking experiment. Scatter plot of smoking *t*-stats in **Figure 3B** can be reproduced with the data contained in this table. 

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Exon_analysis/02_Comparisons.R#L306).



## [Supplementary Table 14](TableS14_DE_jxns_brain_pup_nicotine.tsv) 
**Differentially expressed exon-exon junctions in the nicotine-exposed pup brain.**

Metadata of DE exon-exon junctions in the nicotine pup brain, including for each: 

  * if both the donor and acceptor sites together are known and annotated in GENCODE M25 (`inGencode` variable); 
  * if the donor and acceptor sites are individually annotated in GENCODE M25 (`inGencodeStart` and `inGencodeEnd` variables, respectively); 
  * the junction `Class`: *Novel* (if both start and end sites are unknown, also known as fully novel junctions), *InGen* (already annotated in GENCODE M25), *AltStartEnd* (if it has only one known site), or *ExonSkip* (with sites from non-successive exons, both known individually but not together), and
  * if they are fusion junctions, meaning that they connect exons from different genes (`isFusion` variable). 

  Their logFC, moderated *t*-stats, *p*-value, and adjusted *p*-value are provided. These statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names.

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Jx_analysis/01_Modeling.R#L209).



## [Supplementary Table 15](TableS15_DE_jxns_brain_pup_smoking.tsv)
**Differentially expressed exon-exon junctions in the smoking-exposed pup brain.**

Same as in [**Table S14**](TableS14_DE_jxns_brain_pup_nicotine.tsv) but for DE exon-exon junctions in the smoking pup brain.

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Jx_analysis/01_Modeling.R#L235).



## [Supplementary Table 16](TableS16_DGE_PNE_SVadjusted.tsv)
**Differential gene expression results for nicotine exposure after SV adjustment.**

Gene-level metadata and the logFC, moderated *t*-stats, *p*-value, and adjusted *p*-value for DE of each gene comparing nicotine vs vehicle exposure in the pup brain, after adjusting DGE models for batch-related surrogate variables. The statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names. 

Created [here](https://github.com/LieberInstitute/smoking-nicotine-mouse/blob/9e6e93ce4ea03c1b7aeaa4ea6de73891975422ab/code/04_DEA/sensitivity_analyses/01_litter_and_batch_effects.R#L738).



## [Supplementary Table 17](TableS17_DGE_MSDP_SVadjusted.tsv)
**Differential gene expression results for smoking exposure after SV adjustment.**

Same as in [**Table S16**](TableS16_DGE_PNE_SVadjusted.tsv) but for the smoking exposure.

Created [here](https://github.com/LieberInstitute/smoking-nicotine-mouse/blob/9e6e93ce4ea03c1b7aeaa4ea6de73891975422ab/code/04_DEA/sensitivity_analyses/01_litter_and_batch_effects.R#L741).



## [Supplementary Table 18](TableS18_human_vs_mice_results.tsv)
**Differential gene expression results for gene pairs of mouse-human orthologs.**

The logFC, moderated *t*-stats, *p*-value, and adjusted *p*-value of the human gene for smoking exposure in the prenatal and adult human brain, and for the corresponding mouse orthologous gene for the 5 experimental mice groups (as in [**Table S5**](TableS5_DEA_results_all_genes.tsv)), are presented. Only mouse genes with human ortholog(s) present in the human dataset from [Semick et al., 2020](https://www.nature.com/articles/s41380-018-0223-1) are considered. The DE statistics were computed with `topTable()` from [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html); see its documentation for the definition of the variable names. **Figure 5** and **Figure 18C** were created with the data this table provides. 

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/04_DEA/Gene_analysis/02_Comparisons.R#L937).



## [Supplementary Table 19](TableS19_sample_variable_dict.tsv)
**Dictionary of sample variables.**

Description of the sample variables used throughout this project. 



## [Supplementary Table 20](TableS20_novel_DE_jxns_genes.tsv)
**Associated genes of fully novel DE exon-exon junctions in pup brain.**

Nearest, following, and preceding genes of the fully novel DE exon-exon junctions without assigned gene for the nicotine and smoking exposure in pup brain.

Created [here](https://github.com/LieberInstitute/smokingMouse_Indirects/blob/5bbf5b3e272f7bff90975f23446f281479608818/code/07_Jxn_anno/01_Jxn_anno.R#L309).
