
# Raw data access


## Public access

Through the [_smokingMouse_](https://bioconductor.org/packages/release/data/experiment/html/smokingMouse.html) Bioconductor data package, we offer free open access to all generated and analyzed  [`RangedSummarizedExperiment`](https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html) (RSE) objects containing:

* The raw expression counts and lognorm counts of genes, exons, and exon-exon junctions, and the transcripts per million (TMP) and log-TPM of transcripts, across the 208 mouse samples from pup and adult brain, and adult blood. 
* Sample-level metadata (what's contained in [Maternal_Smoking_pheno.txt](Maternal_Smoking_pheno.txt)) and quality control metrics
* Feature information 

The results of the DEA for smoking exposure in human postmortem prenatal and adult prefrontal cortices from [Semick et al. 2020](https://www.nature.com/articles/s41380-018-0223-1) are also accessible.

See more details of the provided datasets and how to access them in the [_smokingMouse_ ](http://research.libd.org/smokingMouse/) web page.



## Internal access

The raw data from this project lives at different external directories.

The script `code/01_SPEAQeasy/01_make_manifest.R` creates soft links (stored in `raw-data/fastq`) to the original raw data.

Originally, the paths to the original data were determined by reading in: `/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/blood/preprocessed_data_RNAsp/samples.manifest` and `/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/brain/preprocessed_data_RNAsp/samples.manifest`, which ultimately refered to FASTQ files under `/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/FASTQ/`. After data was transferred from dcl01 to dcs05, `code/01_SPEAQeasy/01_make_manifest.R` was re-run, referring to new paths under `/dcs05/lieber/marmaypag/nicotineGonzalez2023_LIBD001`.

SPEAQeasy was run after the first execution of `code/01_SPEAQeasy/01_make_manifest.R` and referenced files under dcl01 (which should not make a difference).
