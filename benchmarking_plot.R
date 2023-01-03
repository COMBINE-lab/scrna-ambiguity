library(ggplot2)
# library(patchwork)

root_dir=file.path("scrna_ambiguity")
f_vec = list.files(root_dir, pattern = "\\.time$", recursive=TRUE)

log_list = sapply(f_vec, function(f) {
    tail(read.table(file.path(root_dir, f), sep = "!"),22)
})

log_df = as.data.frame(do.call(rbind, lapply(log_list, function(x) {
    memory = strsplit(x[9], "\tMaximum resident set size (kbytes): ", fixed = TRUE)[[1]][2]
    memory = as.numeric(memory)/1e6
    wall_clock_time = x[4]
    # cat(wall_clock_time, "\n")
    wall_clock_time = strsplit(wall_clock_time, "\tElapsed (wall clock) time (h:mm:ss or m:ss): ", fixed = TRUE)[[1]][2]
    time_sep_by_colon = strsplit(wall_clock_time, ":")[[1]]
    if (length(time_sep_by_colon) == 2) {
       wall_clock_time = as.numeric(time_sep_by_colon[1]) + as.numeric(time_sep_by_colon[2])/60
    } else if(length(time_sep_by_colon) == 3) {
       wall_clock_time = as.numeric(time_sep_by_colon[2]) + as.numeric(time_sep_by_colon[3])/60 + as.numeric(time_sep_by_colon[1])*60
    } else {
        warning("Needs more conditions")
    }

    c(wall_clock_time, memory)
})))

colnames(log_df) = c("time_minutes", "memory_gbs")

log_df$name = sapply(rownames(log_df), USE.NAMES = FALSE, function(s) {
    filename = basename(file.path(s))
    cmd_name= strsplit(filename, ".time")[[1]][1]
    cmd_name
})

# For the classification experiment, I will compare the time of index building and time of mapping

classification_log_df = log_df[grep("classification_experiment",rownames(log_df)),]

classification_af_quant_log_df = classification_log_df[grepl("piscem_map_sc.time",rownames(classification_log_df)) | 
                                                    grepl("salmon_alevin.time", rownames(classification_log_df)),]

classification_af_quant_log_df$context = sapply(rownames(classification_af_quant_log_df), USE.NAMES = FALSE, function(s) {
    sep_by_slash = strsplit(s, "/", fixed = TRUE)[[1]]
    sep_by_slash[4]
})

classification_af_quant_log_df$method = sapply(rownames(classification_af_quant_log_df), USE.NAMES = FALSE, function(s) {
    sep_by_slash = strsplit(s, "/", fixed = TRUE)[[1]]
    sapply(sep_by_slash[2], function(x) {
        stringr::word(substring(x, first = 12), 1, sep = "_results")
    })
})

classification_af_quant_log_df$name = paste(classification_af_quant_log_df$context, classification_af_quant_log_df$method, sep = "_")


ggplot(data=classification_af_quant_log_df, aes(x=method, y=as.numeric(memory_gbs), fill=context)) +
    geom_bar(position="dodge", stat="identity") + 
  xlab("Method") + ylab("Peak memory in GB") + 
    scale_fill_brewer(palette="Spectral")  
    # + theme(legend.position = "none") 

ggplot(data=classification_af_quant_log_df, aes(x=method, y=as.numeric(time_minutes), fill=context)) +
    geom_bar(position="dodge", stat="identity") + 
  xlab("Method") + ylab("Run time in minutes") + 
    scale_fill_brewer(palette="Spectral")


classification_kb_quant_log_df = classification_log_df[grepl("kallisto_bus.time", rownames(classification_log_df)) ,]


classification_kb_quant_log_df$context = sapply(rownames(classification_kb_quant_log_df), USE.NAMES = FALSE, function(s) {
    sep_by_slash = strsplit(s, "/", fixed = TRUE)[[1]]
    sep_by_slash[4]
})

classification_kb_quant_log_df$method = sapply(rownames(classification_kb_quant_log_df), USE.NAMES = FALSE, function(s) {
    sep_by_slash = strsplit(s, "/", fixed = TRUE)[[1]]
    if(sep_by_slash[5] == "kbd_idx_mature_as_ref_genome_as_dlist") {
        "sc_idx"
    } else {
        "sn_idx"
    }
})

classification_kb_quant_log_df$name = paste(classification_kb_quant_log_df$context, classification_kb_quant_log_df$method, sep = "_")

ggplot(data=classification_kb_quant_log_df, aes(x=method, y=as.numeric(memory_gbs), fill=context)) +
    geom_bar(position="dodge", stat="identity") + 
  xlab("Method") + ylab("Peak memory in GB") + 
    scale_fill_brewer(palette="Spectral")  
    # + theme(legend.position = "none") 

ggplot(data=classification_kb_quant_log_df, aes(x=method, y=as.numeric(time_minutes), fill=context)) +
    geom_bar(position="dodge", stat="identity") + 
  xlab("Method") + ylab("Run time in minutes") + 
    scale_fill_brewer(palette="Spectral")


# for the starsolo simulation
starsim_log_df = log_df[grep("starsolo_simulation",rownames(log_df)),]

starsim_quant_log_df = starsim_log_df[grepl("piscem_map_sc.time",rownames(starsim_log_df)) | 
                                    grepl("kb_count.time", rownames(starsim_log_df)) | 
                                    grepl("simpleaf_quant.time", rownames(starsim_log_df)) | 
                                    grepl("star_quant.time", rownames(starsim_log_df)),]

starsim_quant_log_df[2, "time_minutes"] = sum(starsim_quant_log_df[c(1,2), "time_minutes"])
starsim_quant_log_df[2, "memory_gbs"] = max(starsim_quant_log_df[c(1,2), "memory_gbs"])
starsim_quant_log_df = starsim_quant_log_df[-1, ]

starsim_quant_log_df[grepl("simpleaf_quant.time", rownames(starsim_quant_log_df)), "name"] = stringr::word(substring(basename(dirname(rownames(starsim_quant_log_df)[grepl("simpleaf_quant.time",rownames(starsim_quant_log_df))])),first=4 ), 1, sep = "_quant")
starsim_quant_log_df[!grepl("simpleaf_quant.time", rownames(starsim_quant_log_df)), "name"] = stringr::word(starsim_quant_log_df[!grepl("simpleaf_quant.time", rownames(starsim_quant_log_df)), "name"], 1, sep = "_")

ggplot(data=starsim_quant_log_df, aes(x=name, y=as.numeric(memory_gbs), fill=name)) +
  geom_bar(stat="identity") + 
  xlab("Method") + ylab("Peak memory in GB") + 
    scale_fill_brewer(palette="Spectral")  
    # + theme(legend.position = "none") 

ggplot(data=starsim_quant_log_df, aes(x=name, y=as.numeric(time_minutes), fill=name)) +
  geom_bar(stat="identity") + 
  xlab("Method") + ylab("Run time in minutes") + 
    scale_fill_brewer(palette="Spectral")


# for the mouse brain nuclei datasets
mouse_log_df = log_df[grep("mouse_brain_nuclei",rownames(log_df)),]

mouse_quant_log_df = mouse_log_df[grepl("piscem_map_sc.time",rownames(mouse_log_df)) | 
                                    grepl("kb_count.time", rownames(mouse_log_df)) | 
                                    grepl("simpleaf_quant.time", rownames(mouse_log_df)) | 
                                    grepl("star_quant.time", rownames(mouse_log_df)),]

mouse_adult_quant_log_df = mouse_quant_log_df[grepl("adult_brain",rownames(mouse_quant_log_df)),]

mouse_adult_quant_log_df[2, "time_minutes"] = sum(mouse_adult_quant_log_df[c(1,2), "time_minutes"])
mouse_adult_quant_log_df[2, "memory_gbs"] = max(mouse_adult_quant_log_df[c(1,2), "memory_gbs"])
mouse_adult_quant_log_df = mouse_adult_quant_log_df[-1, ]
mouse_adult_quant_log_df[2, "name"] = "piscem_simpleaf_quant"
mouse_adult_quant_log_df[grepl("simpleaf_quant.time", rownames(mouse_adult_quant_log_df)), "name"] = stringr::word(substring(basename(dirname(dirname(rownames(mouse_adult_quant_log_df)[grepl("simpleaf_quant.time",rownames(mouse_adult_quant_log_df))]))),first=4 ), 1, sep = "_quant")
mouse_adult_quant_log_df[!grepl("simpleaf_quant.time", rownames(mouse_adult_quant_log_df)), "name"] = stringr::word(mouse_adult_quant_log_df[!grepl("simpleaf_quant.time", rownames(mouse_adult_quant_log_df)), "name"], 1, sep = "_")

ggplot(data=mouse_adult_quant_log_df, aes(x=name, y=as.numeric(memory_gbs), fill=name)) +
  geom_bar(stat="identity") + 
  xlab("Method") + ylab("Peak memory in GB") + 
    scale_fill_brewer(palette="Spectral")  
    # + theme(legend.position = "none") 

ggplot(data=mouse_adult_quant_log_df, aes(x=name, y=as.numeric(time_minutes), fill=name)) +
  geom_bar(stat="identity") + 
  xlab("Method") + ylab("Run time in minutes") + 
    scale_fill_brewer(palette="Spectral")



mouse_E18_quant_log_df = mouse_quant_log_df[grepl("chsz",rownames(mouse_quant_log_df)),]

mouse_E18_quant_log_df[2, "time_minutes"] = sum(mouse_E18_quant_log_df[c(1,2), "time_minutes"])
mouse_E18_quant_log_df[2, "memory_gbs"] = max(mouse_E18_quant_log_df[c(1,2), "memory_gbs"])
mouse_E18_quant_log_df = mouse_E18_quant_log_df[-1, ]
mouse_E18_quant_log_df[2, "name"] = "piscem_simpleaf_quant"
mouse_E18_quant_log_df[grepl("simpleaf_quant.time", rownames(mouse_E18_quant_log_df)), "name"] = stringr::word(substring(basename(dirname(dirname(rownames(mouse_E18_quant_log_df)[grepl("simpleaf_quant.time",rownames(mouse_E18_quant_log_df))]))),first=4 ), 1, sep = "_quant")
mouse_E18_quant_log_df[!grepl("simpleaf_quant.time", rownames(mouse_E18_quant_log_df)), "name"] = stringr::word(mouse_E18_quant_log_df[!grepl("simpleaf_quant.time", rownames(mouse_E18_quant_log_df)), "name"], 1, sep = "_")

ggplot(data=mouse_E18_quant_log_df, aes(x=name, y=as.numeric(memory_gbs), fill=name)) +
  geom_bar(stat="identity") + 
  xlab("Method") + ylab("Peak memory in GB") + 
    scale_fill_brewer(palette="Spectral")  
    # + theme(legend.position = "none") 

ggplot(data=mouse_E18_quant_log_df, aes(x=name, y=as.numeric(time_minutes), fill=name)) +
  geom_bar(stat="identity") + 
  xlab("Method") + ylab("Run time in minutes") + 
    scale_fill_brewer(palette="Spectral")




