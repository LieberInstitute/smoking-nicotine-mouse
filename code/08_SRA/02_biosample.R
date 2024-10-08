library(tidyverse)
library(here)
library(sessioninfo)

man_path = here('processed-data', '01_SPEAQeasy', 'samples.manifest')
pheno_path = here('raw-data', 'Maternal_Smoking_pheno.txt')
out_path = here('processed-data', '08_SRA', 'biosample.tsv')

#   Read in the manifest just to make sure the phenotype table includes the
#   expected set of IDs (though we'll only really use the phenotype table)
meta_df = read.table(
        man_path, col.names = c('r1', 'md1', 'r2', 'md2', 'sample_name')
    ) |>
    as_tibble() |>
    select(sample_name)

pheno_df = read_tsv(pheno_path, show_col_types = FALSE) |>
    #   Match column names expected by SRA
    rename(dev_stage = Age, tissue = Tissue) |>
    mutate(
        #   Some sample IDs differ in name for the same samples between
        #   'pheno_df' and 'meta_df'
        sample_name = str_replace(
            SAMPLE_ID, '^(F.*?)-blood$', 'Sample_\\1_blood'
        ),
        organism = 'mus musculus',
        isolate = sprintf('%s-exposed %s group', Expt, tolower(Group)),
        collection_date = ifelse(
            is.na(date),
            'not collected',
            sprintf(
                '20%s-%s-%s',
                substr(date, 5, 6), substr(date, 1, 2), substr(date, 3, 4)
            )
        ),
        geo_loc_name = 'United States: Baltimore, MD',
        sex = ifelse(Sex == "M", 'male', 'female'),
        treatment = Expt
    ) |>
    select(
        sample_name, organism, isolate, dev_stage, collection_date,
        geo_loc_name, sex, tissue, treatment
    )

#   Phenotype data and manifest should include same set of sample IDs
stopifnot(setequal(pheno_df$sample_name, meta_df$sample_name))

write_tsv(pheno_df, out_path)

session_info()
