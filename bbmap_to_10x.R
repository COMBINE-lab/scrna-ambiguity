#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

sc_m_r2_path = args[1]
sc_n_r2_path = args[2]
sn_m_r2_path = args[3]
sn_n_r2_path = args[4]

sc_m_r1_path = args[5]
sc_n_r1_path = args[6]
sn_m_r1_path = args[7]
sn_n_r1_path = args[8]
num_cores = as.integer(args[9])

suppressPackageStartupMessages({
    library(Biostrings)
    library(doParallel)
})

cat("Start generating random read1\n")

# generate random CB+UMI

registerDoParallel(cores=num_cores)
set.seed = 17
r1_size = 28
n_barcode = 10000000
batch_size= ceiling(n_barcode/num_cores)
rand_r1 <- DNAStringSet(unique(foreach(i=1:num_cores, .combine=c) %dopar% {
    ab = c("A", "C", "T", "G")
    sapply(1:batch_size, function(x) paste(sample(ab,  size = r1_size, replace = TRUE),  collapse = ""))
}))

registerDoSEQ()

cat("Start reading read2 files\n")

# read the simulated single-end reads
sc_m_r2 = readDNAStringSet(sc_m_r2_path, format = "fastq")
sc_n_r2 = readDNAStringSet(sc_n_r2_path, format = "fastq")
sn_m_r2 = readDNAStringSet(sn_m_r2_path, format = "fastq")
sn_n_r2 = readDNAStringSet(sn_n_r2_path, format = "fastq")

cat("Start assigning names to read1\n")

# Assign read 1 to each read
## sc mature
first_id = 1
last_id = length(sc_m_r2)

sc_m_r1 = rand_r1[first_id:last_id]
names(sc_m_r1) = names(sc_m_r2)

## sc nascent
first_id = last_id + 1
last_id = length(sc_m_r2) + length(sc_n_r2)


sc_n_r1 = rand_r1[first_id:last_id]
names(sc_n_r1) = names(sc_n_r2)

## sn mature
first_id = 1
last_id = length(sn_m_r2)
sn_m_r1 = rand_r1[first_id:last_id]
names(sn_m_r1) = names(sn_m_r2)

## sn nascent
first_id = last_id + 1
last_id = length(sn_m_r2) + length(sn_n_r2)

sn_n_r1 = rand_r1[first_id:last_id]
names(sn_n_r1) = names(sn_n_r2)


cat("Start writing read1 files\n")

# write files
# read1
# sc mature - 4000000
writeXStringSet(sc_m_r1, sc_m_r1_path, compress = FALSE, format="fastq")

# sc nascent - 1000000
writeXStringSet(sc_n_r1, sc_n_r1_path, compress = FALSE, format="fastq")

# sn mature - 1000000
writeXStringSet(sn_m_r1, sn_m_r1_path, compress = FALSE, format="fastq")

# sn nascent - 4000000
writeXStringSet(sn_n_r1, sn_n_r1_path, compress = FALSE, format="fastq")
