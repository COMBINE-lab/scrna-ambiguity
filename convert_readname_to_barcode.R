args = commandArgs(trailingOnly=TRUE)
# args = c()
# args[1]="/mnt/scratch2/scrna_ambiguity/classification_experiment/data/bbmap_sim/matching_based_truth"
# args[2]="/mnt/scratch2/scrna_ambiguity/classification_experiment/data/bbmap_sim/read_fastqs"
suppressPackageStartupMessages({
    library(Biostrings)
    library(GenomicFeatures)
})
matching_based_truth_dir = args[1]
read_fastqs_dir = args[2]

simulation_names = c("cytoplasmic_mature", "cytoplasmic_nascent", "nuclear_mature", "nuclear_nascent")


get_name_to_bc <- function(read1_path) {
    read1 = readDNAStringSet(read1_path, format = "fastq")
    name_to_bc = as.character(read1)

    # the empirical splice status program keeps only the read name before the space in it.
    names(name_to_bc) = stringr::word(names(read1), 1)
    name_to_bc
}

for (s_name in simulation_names) {
    # get read name to bc mapping
    read1_path = file.path(read_fastqs_dir, s_name, paste0(s_name, "_R1.fastq"))
    name_to_bc = get_name_to_bc(read1_path)

    # get matching based truth
    matching_based_truth_path = file.path(matching_based_truth_dir, paste0(s_name, "_splice_status_truth.tsv"))
    matching_based_truth = read.csv(matching_based_truth_path,
                                    header=FALSE, sep="\t"
                                    )

    # convert read name to barcode and write down
    matching_based_truth$V1 = name_to_bc[matching_based_truth$V1]
    write.table(matching_based_truth, 
                matching_based_truth_path, 
                col.names=FALSE, 
                row.names=FALSE, 
                quote=FALSE,
                sep="\t")
}
