library(tidyverse)
library(here)
library(sessioninfo)

man_path = here('processed-data', '01_SPEAQeasy', 'samples.manifest')

meta_df = read.table(
        man_path, col.names = c('r1', 'md1', 'r2', 'md2', 'sample_name')
    ) |>
    as_tibble() |>
    select(c('filename', 'filename2', 'sample_name'))

#   Show that a sample always consists of exactly 2 files
total_files = meta_df |>
    pivot_longer(cols = c('filename', 'filename2')) |>
    group_by(sample_id) |>
    summarize(num_files = n()) |>
    pull(num_files)
stopifnot(all(total_files == 2))

#   Re-order columns to match SRA's expectations
meta_df = meta_df |>
    select(
        c(
            'sample_name', 'library_ID', 'title', 'library_strategy',
            'library_source', 'library_selection', 'library_layout', 'platform',
            'instrument_model', 'design_description', 'filetype', 'filename',
            'filename2'
        )
    )

session_info()
