#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

genome_path = args[1]
gene_gtf_path = args[2]
spliceu_dir_path = args[3]
filename_prefix = args[4]

if (is.na(filename_prefix)) {
    filename_prefix = "spliceu"
}

suppressPackageStartupMessages({
    library(Biostrings)
    library(GenomicFeatures)
    library(BSgenome)
})

genome = readDNAStringSet(genome_path, format = "fasta")
names(genome) = stringr::word(names(genome), 1)

# mature transcripts
suppressWarnings({
    suppressMessages({
        txdb <- GenomicFeatures::makeTxDbFromGFF(gene_gtf_path, format = "gtf")
        mature_txps = extractTranscriptSeqs(genome, txdb, use.names=TRUE)
        writeXStringSet(mature_txps, file.path(spliceu_dir_path, paste0(filename_prefix ,"_mature.fa")))

        t2g_mature <- AnnotationDbi::select(txdb, keys = transcripts(txdb)$tx_name, 
                                        keytype = "TXNAME", columns = "GENEID")
        t2g_mature$STATUS = "S"
    })
}) 


# nascent transcripts
nascent_txps = getSeq(genome, genes(txdb))
writeXStringSet(nascent_txps, file.path(spliceu_dir_path, paste0(filename_prefix ,"_nascent.fa")))

# gene id doesn't have -I, txp id has -I
t2g_nascent = data.frame(TXNAME = paste0(names(nascent_txps), "-I"), GENEID = names(nascent_txps), STATUS = "U")
names(nascent_txps) = t2g_nascent$TXNAME

# write 
dir.create(spliceu_dir_path, showWarnings=FALSE, recursive=TRUE)

# write sequences
writeXStringSet(c(mature_txps, nascent_txps), file.path(spliceu_dir_path, paste0(filename_prefix ,"_all.fa")))

# write t2g after combining mature and nascent
write.table(rbind(t2g_mature, t2g_nascent), file.path(spliceu_dir_path, paste0(filename_prefix ,"_t2g_3col.tsv")), col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")
write.table(rbind(t2g_mature[,c("TXNAME","GENEID")], t2g_nascent[,c("TXNAME","GENEID")]), file.path(spliceu_dir_path, paste0(filename_prefix ,"_t2g.tsv")), col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")
write.table(unique(cbind(t2g_nascent[,"GENEID"], t2g_nascent[,"GENEID"])), file.path(spliceu_dir_path, paste0(filename_prefix ,"_g2g.tsv")), col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")

