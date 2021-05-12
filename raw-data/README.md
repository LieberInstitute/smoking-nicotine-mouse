The raw data from this project lives at different external directories.

The script `code/01_SPEAQeasy/01_make_manifest.R` creates soft links to the original raw data. The paths to the original data are determined by reading in `/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/blood/preprocessed_data_RNAsp/samples.manifest` and `/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/brain/preprocessed_data_RNAsp/samples.manifest`, which ultimately refer to FASTQ files under `/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/FASTQ/`.

SPEAQeasy then references the soft links in `raw-data/fastq` to determine the inputs to the pipeline.
