library('here')
library('jaffelab')
library('SummarizedExperiment')

man_path = here("processed-data/01_SPEAQeasy/samples.manifest")
rse_dir = here("processed-data/01_SPEAQeasy/pipeline_output/count_objects")

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

#   Load RSE objects and add flowcell as a column to the colData
rse_gene_path = here(rse_dir, "rse_gene_smoking_mouse_n208.Rdata")
rse_exon_path = here(rse_dir, "rse_exon_smoking_mouse_n208.Rdata")
rse_jx_path = here(rse_dir, "rse_jx_smoking_mouse_n208.Rdata")
rse_tx_path = here(rse_dir, "rse_tx_smoking_mouse_n208.Rdata")

load(rse_gene_path)
colData(rse_gene)$flowcell = flowcells

load(rse_exon_path)
colData(rse_exon)$flowcell = flowcells

load(rse_jx_path)
colData(rse_jx)$flowcell = flowcells

load(rse_tx_path)
colData(rse_tx)$flowcell = flowcells
