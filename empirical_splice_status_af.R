#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


get_read_splice_status_af <- function(map_txt_path) {
    ## parse mapping records
    read_map = read.csv(map_txt_path, header=FALSE, sep = "\t")
    colnames(read_map) = c("ID", "HI", "NH", "CB", "UMI", "DIR", "REF_I")
    read_map$CBUMI = paste0(substring(read_map$CB, first=4), substring(read_map$UMI, first=5))
    read_map = read_map[,-c(1,2,4:6)]
    read_map$NH = as.integer(substring(read_map$NH, first=4))
    read_map$REF_I = grepl("-I", read_map$REF_I)

    # single-cell mature
    # aggregate REF_I (if nascent mapping) by CBUMI 
    read_map_agg = aggregate(REF_I~., data = read_map, sum)

    # assign splice status to each read 
    ## init vector by ambiguous
    read_splice_status_af = rep("a", nrow(read_map_agg))
    names(read_splice_status_af) = read_map_agg$CBUMI

    ## if no nascent mapping, it's mature
    read_splice_status_af[read_map_agg$REF_I == 0] = "m"

    ## if all nascent mappings, it's nascent
    read_splice_status_af[read_map_agg$REF_I == read_map_agg$NH] = "n"

    read_splice_status_af
}

map_txt_path = args[1]
out_path = args[2]

read_splice_status = get_read_splice_status_af(map_txt_path)

write.table(cbind(names(read_splice_status), read_splice_status), out_path, col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")



