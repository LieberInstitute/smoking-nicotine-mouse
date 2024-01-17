The raw data from this project lives at different external directories.

The script `code/01_SPEAQeasy/01_make_manifest.R` creates soft links (stored in `raw-data/fastq`) to the original raw data. Originally, the paths to the original data were determined by reading in `/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/blood/preprocessed_data_RNAsp/samples.manifest` and `/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/brain/preprocessed_data_RNAsp/samples.manifest`, which ultimately refered to FASTQ files under `/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/FASTQ/`. After data was transferred from dcl01 to dcs05, `code/01_SPEAQeasy/01_make_manifest.R` was re-run, referring to new paths under `/dcs05/lieber/marmaypag/nicotineGonzalez2023_LIBD001`.

SPEAQeasy was run after the first execution of `code/01_SPEAQeasy/01_make_manifest.R` and referenced files under dcl01 (which should not make a difference).
