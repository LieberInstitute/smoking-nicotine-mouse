library(tidyverse)
library(here)
library(sessioninfo)

man_path = here('processed-data', '01_SPEAQeasy', 'samples.manifest')
pheno_path = here('raw-data', 'Maternal_Smoking_pheno.txt')
out_path = here('processed-data', '08_SRA', 'metadata.tsv')
fastq_dir = here('raw-data', 'fastq')
design_description_text = paste(
    "Animals",
    "Wild type male and female mice (C57BL/6J; stock # 000664, Jackson Laboratories) were purchased and used for timed breeding. 6 week old female mice were paired with male mice. Copulation plugs were checked daily and male mice were removed upon identification of plugs. Female mice were monitored for pregnancy, and separated upon pregnancy confirmation. Pups were euthanized by decapitation on the first day following birth, e.g. postnatal day 0 (P0).  Pregnant dams were euthanized by decapitation and trunk blood was collected into a heparinized tube. Brains were rapidly extracted from the skull and the frontal cortex was dissected from the brain over wet ice on a steel block using a scalpel. Frontal cortical tissue was snap frozen in chilled 2-methylbutane. Samples were transferred to tubes and placed on dry ice and stored at -80°C until further processing for RNA extraction. All experiments and procedures were approved by the Johns Hopkins Animal Care and Use Committee and in accordance with the Guide for the Care and Use of Laboratory Animals.",
    "Nicotine administration",
    "Free-base(-)-nicotine (Sigma) was dissolved in normal saline. Nicotine (1.5mg/kg) or vehicle (saline) was administered to female dams (2X/daily - 8AM and 4PM). Administration started the week before mice were paired and continued until E17.",
    "Smoking exposure",
    "Pregnant dams were placed into a smoking chamber for 5 hours/day, 5 days/week starting one week before mice were paired for breeding and until the time of delivery. This chamber contains a smoking machine (Model TE-10, Teague Enterprises, Davis, CA) that burns 5 cigarettes (2R4F reference cigarettes (2.45 mg nicotine/cigarette; Tobacco Research Institute, University of Ky) at a time, taking 2 second duration puffs at a flow rate of 1.05 l/min, to provide a standard puff of 35 cm3, providing a total of 8 puffs per minute. The machine is adjusted to produce side stream (89%) and mainstream smoke (11%). The chamber atmosphere is monitored to maintain total suspended particulate at 90 mg/m3, and carbon monoxide at 350 ppm. Control pregnant dams were kept in a filtered air environment.",
    "Tissue processing and RNA isolation and sequencing",
    "Total RNA was extracted from samples using Trizol followed by purification with an RNeasy Micro kit (Qiagen). Paired-end strand-specific sequencing libraries were prepared and sequenced by Macrogen from 1ug total RNA using the TruSeq Stranded mRNA kit with ERCC Spike in. Libraries were sequenced on an Illumina NovaSeq6000 S4, 150bp paired end. Output was targeted at 60M total reads (R1 30M and R2 30M) million 150-bp paired-end reads."
)

meta_df = read.table(
        man_path,
        col.names = c('filename', 'md1', 'filename2', 'md2', 'sample_name')
    ) |>
    as_tibble() |>
    select(filename, filename2, sample_name) |>
    left_join(
        read_tsv(pheno_path, show_col_types = FALSE) |>
            mutate(
                #   Some sample IDs differ in name for the same samples between
                #   pheno data and meta_df
                sample_name = str_replace(
                    SAMPLE_ID, '^(F.*?)-blood$', 'Sample_\\1_blood'
                )
            ),
        by = 'sample_name'
    ) |>
    mutate(
        library_ID = sample_name,
        title = sprintf(
            'Bulk RNA-seq of mouse %s: %s-exposure %s group',
            tolower(Tissue), tolower(Expt), tolower(Group)
        ),
        library_strategy = 'RNA-Seq',
        library_source = 'TRANSCRIPTOMIC',
        library_selection = 'PolyA',
        library_layout = 'paired',
        platform = 'ILLUMINA',
        instrument_model = 'Illumina NovaSeq 6000',
        design_description = factor(design_description_text),
        filetype = 'fastq'
    ) |>
    #   Re-order columns to match SRA's expectations
    select(
        sample_name, library_ID, title, library_strategy,
        library_source, library_selection, library_layout, platform,
        instrument_model, design_description, filetype, filename,
        filename2
    )

#   We actually want to use the FASTQ soft links under 'fastq_dir' for upload
#   since the directory structure is flat (an SRA requirement). Form a table
#   of links and destination files
fastq_mapping = tibble(
        link_path = list.files(fastq_dir, full.names = TRUE)
    ) |>
    mutate(dest_path = Sys.readlink(link_path))

#   'meta_df' (from the SPEAQeasy manifest) should refer to the same set of
#   FASTQs as are present under 'fastq_dir'
stopifnot(
    setequal(fastq_mapping$dest_path, c(meta_df$filename, meta_df$filename2))
)

#   Use the symbolic links, not actual files
meta_df = meta_df |>
    mutate(
        filename = fastq_mapping$link_path[
            match(filename, fastq_mapping$dest_path)
        ],
        filename2 = fastq_mapping$link_path[
            match(filename2, fastq_mapping$dest_path)
        ]
    )

#   Also check disk usage (to make sure we comply with SRA requirements)
disk_usage = sapply(
    c(meta_df$filename, meta_df$filename2),
    function(x) {
        as.integer(
            system(sprintf('du -k $(readlink %s) | cut -f 1', x), intern = TRUE)
        )
    }
)
message(
    sprintf(
        "FASTQ files occupy a total of %sGB.", round(sum(disk_usage) / 1e6, 1)
    )
)

#   Individual files must be less than 100GB
stopifnot(all(disk_usage < 100e6))

write_tsv(meta_df, out_path)

session_info()
