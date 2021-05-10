#  This script was run interactively to generate 'samples.manifest' and link
#  existing FASTQ files into the 'raw-data/fastq' directory.

library('here')

#  Paths, including hard paths to the original manifests
blood_man_path = '/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/blood/preprocessed_data_RNAsp/samples.manifest'
brain_man_path = '/dcl01/lieber/ajaffe/lab/Nicotine/nicotine_mouse_rnaseq/brain/preprocessed_data_RNAsp/samples.manifest'
out_dir = file.path(here::here(), 'raw-data', 'fastq')
new_man_path = file.path(here::here(), 'processed-data', '01_SPEAQeasy', 'samples.manifest')

#  Read in original manifests into one table, and verify files exist
orig_man = rbind(read.table(blood_man_path), read.table(brain_man_path))
stopifnot(all(file.exists(orig_man[,1])))
stopifnot(all(file.exists(orig_man[,3])))

#  Create links in this repository to the original FASTQ files
dir.create(out_dir, showWarnings=FALSE)
file.symlink(orig_man[,1], file.path(out_dir, basename(orig_man[,1])))
file.symlink(orig_man[,3], file.path(out_dir, basename(orig_man[,3])))

#  Create a new manifest using the paths to the FASTQ links
new_man = paste(file.path(out_dir, basename(orig_man[,1])),
                0,
                file.path(out_dir, basename(orig_man[,3])),
                0,
                orig_man[,5],
                sep='\t')

dir.create(dirname(new_man_path), showWarnings=FALSE)
writeLines(new_man, con=new_man_path)
