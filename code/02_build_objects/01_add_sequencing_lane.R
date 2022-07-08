library('here')
library('jaffelab')
library('SummarizedExperiment')

man_path = here("processed-data/01_SPEAQeasy/samples.manifest")
out_path = here("processed-data", "build_objects", "flowcell_info.tsv")

#   Read in the manifest to extract the FASTQ paths for the first read in each
#   sample
man = read.table(man_path)
fq_r1 = man[,1]

#   For the first read of each sample, grab the first sequence header, which
#   contains information about the flowcell for each sample
headers = sapply(
    fq_r1,
    function(x) {
        system(paste('gunzip -c', x, '| grep "^@" | head -n 1'), intern = TRUE)
    }
)

#   Grab only the flowcell, discarding other information
flowcells = ss(headers, ':', 3)

#   There are 6 unique flowcells used in the experiment
print('Unique flowcell IDs used in this experiment:')
print(unique(flowcells))

#   Write a table containing sample IDs and associated flowcells
flowcell_table = data.frame('sample_id' = man[,5], 'flowcell' = flowcells)
write.table(
    flowcell_table, file = out_path, quote = FALSE, sep = '\t',
    row.names = FALSE
)

#   This table can then be read in with 'read.table(out_path, header = TRUE)'
