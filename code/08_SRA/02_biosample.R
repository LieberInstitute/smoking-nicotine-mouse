library(tidyverse)
library(here)
library(sessioninfo)

man_path = here('processed-data', '01_SPEAQeasy', 'samples.manifest')
pheno_path = here('raw-data', 'Maternal_Smoking_pheno.txt')

#   Read in the manifest just to make sure the phenotype table includes the
#   expected set of IDs (though we'll only really use the phenotype table)
meta_df = read.table(
        man_path, col.names = c('r1', 'md1', 'r2', 'md2', 'sample_name')
    ) |>
    as_tibble() |>
    select(sample_name)

pheno_df = read.table(pheno_path, header = TRUE) |>
    as_tibble() |>
    #   Match column names expected by SRA
    rename(
        `Sample Name` = SAMPLE_ID,
        age = Age,
        `collection date` = date,
        sex = Sex,
        tissue = Tissue
    ) |>
    mutate(
        Organism = 'mus musculus',
        `geographic location` = 'Baltimore, MD, USA',
        #   Some sample IDs differ in name for the same samples between
        #   'pheno_df' and 'meta_df'
        `Sample Name` = str_replace(
            `Sample Name`, '^(F.*?)-blood$', 'Sample_\\1_blood'
        )
    )

#   Phenotype data and manifest should include same set of sample IDs
stopifnot(
    identical(
        sort(unique(pheno_df[['Sample Name']])),
        sort(meta_df$sample_name)
    )
)

session_info()
