#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

get_read_splice_status_kb <- function(out_bus_sc_index_path, out_bus_sn_index_path) {
    ## bus file
    ## sc index bus
    out_bus_sc_index = read.csv(out_bus_sc_index_path, header=FALSE, sep = "\t")
    colnames(out_bus_sc_index) = c("CB", "UMI", "ECID", "COUNT")
    sc_index_mapped_reads = paste0(out_bus_sc_index$CB, out_bus_sc_index$UMI)

    ## sn index bus
    out_bus_sn_index = read.csv(out_bus_sn_index_path, header=FALSE, sep = "\t")
    colnames(out_bus_sn_index) = c("CB", "UMI", "ECID", "COUNT")
    sn_index_mapped_reads = paste0(out_bus_sn_index$CB, out_bus_sn_index$UMI)

    # assign splice status to each read
    mapped_reads = union(sc_index_mapped_reads, sn_index_mapped_reads)
    read_splice_status = rep("a", length(mapped_reads))
    names(read_splice_status) = mapped_reads
    read_splice_status[setdiff(sc_index_mapped_reads, sn_index_mapped_reads)] = "m"
    read_splice_status[setdiff(sn_index_mapped_reads, sc_index_mapped_reads)] = "n"
    
    read_splice_status
}

out_bus_sc_index_path = args[1]
out_bus_sn_index_path = args[2]
out_path = args[3]

read_splice_status = get_read_splice_status_kb(out_bus_sc_index_path, out_bus_sn_index_path)

write.table(cbind(names(read_splice_status), read_splice_status), out_path, col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")



