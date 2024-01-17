#  This script was run interactively to generate 'samples.manifest' and link
#  existing FASTQ files into the 'raw-data/fastq' directory.
#
#  NOTE: due to a poorly specified gitignore, the original manifest was not
#  committed and therefore lost. The original paths on dcl01 became
#  inaccessible, breaking symlinks, and thus made it necessary to re-link
#  the files. This script was re-run with the "fixed" (dcl05) paths, and thus
#  the original manifest (used in SPEAQeasy) differs from the final version
#  produced on the second run of this script. 

library('here')
library('tidyverse')

#   Because of the transfer from dcl01 to dcs05
old_nicotine_dir = '/dcl01/lieber/ajaffe/lab/Nicotine'
new_nicotine_dir = '/dcs05/lieber/marmaypag/nicotineGonzalez2023_LIBD001'

#  Paths, including hard paths to the original manifests (originally starting
#  with [old_nicotine_dir] instead of [new_nicotine_dir]!)
blood_man_path = file.path(
    new_nicotine_dir,
    'nicotine_mouse_rnaseq/blood/preprocessed_data_RNAsp/samples.manifest'
)
brain_man_path = file.path(
    new_nicotine_dir,
    'nicotine_mouse_rnaseq/brain/preprocessed_data_RNAsp/samples.manifest'
)
out_dir = here('raw-data', 'fastq')
new_man_path = here('processed-data', '01_SPEAQeasy', 'samples.manifest')

#  Read in original manifests into one table, and verify files exist
orig_man = rbind(read.table(blood_man_path), read.table(brain_man_path)) |>
    as_tibble()
colnames(orig_man) = c('r1', 'md1', 'r2', 'md2', 'sample_id')
orig_man = orig_man |>
    mutate(
        across(
            matches('^r[12]$'),
            ~ str_replace(.x, old_nicotine_dir, new_nicotine_dir)
        )
    )
stopifnot(all(file.exists(orig_man$r1)))
stopifnot(all(file.exists(orig_man$r2)))

#  Create links in this repository to the original FASTQ files
dir.create(out_dir, showWarnings=FALSE)
all(file.symlink(orig_man$r1, file.path(out_dir, basename(orig_man$r1))))
all(file.symlink(orig_man$r2, file.path(out_dir, basename(orig_man$r2))))

#  Write a new manifest using the paths to the FASTQ links
dir.create(dirname(new_man_path), showWarnings=FALSE)
write.table(
    orig_man, file = new_man_path, sep = '\t', quote = FALSE, row.names = FALSE,
    col.names = FALSE
)
