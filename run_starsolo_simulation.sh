#!/usr/bin/env bash

script_dir="$1"

. $script_dir/run_me.config

# make working dir
starsim_data_dir="$starsim_root_dir/data"
mkdir -p $starsim_data_dir
#---------------------------------------------------------------------------------------------------------------#
echo "  - Preparing STARsolo simulation reads"
# prepare data
starsim_read_dir="$starsim_data_dir/starsim_reads"

## simulation reads
cmd="wget -qO- https://umd.box.com/shared/static/d2hzsqtqqqn285xal8lmqe7dy2va2jir.gz | tar xzf - -C $starsim_data_dir"
eval $cmd

starsim_read1_path="$starsim_read_dir/_R1_.fq"
starsim_read2_path="$starsim_read_dir/_R2_.fq"
starsim_truth_dir="$starsim_dir/truth"

# reference
echo "  - Fetching human CR3 reference"
cmd="wget -qO- https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz | tar xzf - -C $starsim_data_dir"
eval $cmd

humanCR3_ref_dir="$starsim_data_dir/refdata-cellranger-GRCh38-3.0.0"
humanCR3_genome_path="$humanCR3_ref_dir/fasta/genome.fa"
humanCR3_genes_path="$humanCR3_ref_dir/genes/genes.gtf"

# make spliceu reference
echo "  - Making spliceu reference"
humanCR3_spliceu_ref_dir="$humanCR3_ref_dir/refdata-cellranger-GRCh38-3.0.0_spliceu_ref"
humanCR3_spliceu_ref_prefix="refdata-cellranger-GRCh38-3.0.0_spliceu"
mkdir -p $humanCR3_spliceu_ref_dir

cmd="Rscript $script_dir/make_spliceu_txome.R $humanCR3_genome_path $humanCR3_genes_path $humanCR3_spliceu_ref_dir $humanCR3_spliceu_ref_prefix"
eval $cmd

humanCR3_spliceu_mature_path="$humanCR3_spliceu_ref_dir/${humanCR3_spliceu_ref_prefix}_mature.fa"
humanCR3_spliceu_nascent_path="$humanCR3_spliceu_ref_dir/${humanCR3_spliceu_ref_prefix}_nascent.fa"
humanCR3_spliceu_all_path="$humanCR3_spliceu_ref_dir/${humanCR3_spliceu_ref_prefix}_all.fa"
humanCR3_spliceu_t2g_3col_path="$humanCR3_spliceu_ref_dir/${humanCR3_spliceu_ref_prefix}_t2g_3col.tsv"

# 2. Quantification 
#---------------------------------------------------------------------------------------------------------------#
echo "  - Running alevin-fry splici"
## alevin-fry splici
starsim_af_splici_dir="$starsim_root_dir/alevin_fry_splici_results"
mkdir -p $starsim_af_splici_dir

### splici
starsim_af_splici_idx_dir="$starsim_af_splici_dir/af_splici_idx"
mkdir -p $starsim_af_splici_idx_dir

# build index
cmd="$time -v $simpleaf index -o $starsim_af_splici_idx_dir -t $n_threads -f $humanCR3_genome_path -g $humanCR3_genes_path -r 91 > $starsim_af_splici_idx_dir/simpleaf_index.time 2>&1"
eval $cmd

# quant
starsim_af_splici_quant_dir="$starsim_af_splici_dir/af_splici_quant"
mkdir -p $starsim_af_splici_quant_dir
cmd="$time -v $simpleaf quant -c 10xv3 -o $starsim_af_splici_quant_dir -t $n_threads -i $starsim_af_splici_idx_dir/index  -u -r cr-like -m $starsim_af_splici_idx_dir/index/t2g_3col.tsv -1 $starsim_read1_path -2 $starsim_read2_path > $starsim_af_splici_quant_dir/simpleaf_quant.time 2>&1"
eval $cmd


#---------------------------------------------------------------------------------------------------------------#
echo "  - Running alevin-fry spliceu"
### spliceu
starsim_af_spliceu_dir="$starsim_root_dir/alevin_fry_spliceu_results"
mkdir -p $starsim_af_spliceu_dir

starsim_af_spliceu_idx_dir="$starsim_af_spliceu_dir/af_spliceu_idx"
mkdir -p $starsim_af_spliceu_idx_dir

# build index
cmd="$time -v $simpleaf index --keep-duplicates -o $starsim_af_spliceu_idx_dir -t $n_threads --ref-seq $humanCR3_spliceu_all_path > $starsim_af_spliceu_idx_dir/simpleaf_index.time 2>&1"
eval $cmd

# quant
starsim_af_spliceu_quant_dir="$starsim_af_spliceu_dir/af_spliceu_quant"
mkdir -p $starsim_af_spliceu_quant_dir
cmd="$time -v $simpleaf quant -c 10xv3 -o $starsim_af_spliceu_quant_dir -t $n_threads -i $starsim_af_spliceu_idx_dir/index  -u -r cr-like -m $humanCR3_spliceu_t2g_3col_path -1 $starsim_read1_path -2 $starsim_read2_path > $starsim_af_spliceu_quant_dir/simpleaf_quant.time 2>&1"
eval $cmd

#---------------------------------------------------------------------------------------------------------------#
echo "  - Running alevin-fry spliceu piscem"
# spliceu piscem
starsim_af_spliceu_piscem_dir="$starsim_root_dir/alevin_fry_spliceu_piscem_results"
mkdir -p $starsim_af_spliceu_piscem_dir

# build index
starsim_af_spliceu_piscem_idx_dir="$starsim_af_spliceu_piscem_dir/af_spliceu_piscem_idx"
mkdir -p $starsim_af_spliceu_piscem_idx_dir

cmd="$time -v $piscem build -s $humanCR3_spliceu_all_path -k 31 -m 19 -t $n_threads -o $starsim_af_spliceu_piscem_idx_dir/af_spliceu_piscem > $starsim_af_spliceu_piscem_idx_dir/piscem_build.time 2>&1"
eval $cmd

# quant
## adult brain
starsim_af_spliceu_piscem_quant_dir="$starsim_af_spliceu_piscem_dir/af_spliceu_piscem_quant"
mkdir -p $starsim_af_spliceu_piscem_quant_dir

# mapping
starsim_af_spliceu_piscem_map_dir="$starsim_af_spliceu_piscem_quant_dir/piscem_map"
mkdir -p $starsim_af_spliceu_piscem_map_dir

cmd="$time -v $piscem map-sc -i $starsim_af_spliceu_piscem_idx_dir/af_spliceu_piscem -g chromium_v3 -1 $starsim_read1_path -2 $starsim_read2_path -t $n_threads -o $starsim_af_spliceu_piscem_map_dir > $starsim_af_spliceu_piscem_map_dir/piscem_map_sc.time 2>&1"
eval $cmd

# quantification
cmd="$time -v $simpleaf quant -c 10xv3 -o $starsim_af_spliceu_piscem_quant_dir -t $n_threads --map-dir $starsim_af_spliceu_piscem_map_dir -u -r cr-like -m $humanCR3_spliceu_t2g_3col_path > $starsim_af_spliceu_piscem_quant_dir/simpleaf_quant.time 2>&1"
eval $cmd

#---------------------------------------------------------------------------------------------------------------#
echo "  - Running STARsolo"
# STARsolo
starsim_star_dir="$starsim_root_dir/starsolo_results"
mkdir -p $starsim_star_dir

# build index
cd $starsim_star_dir
starsim_star_idx_dir="$starsim_star_dir/star_index"
mkdir -p $starsim_star_idx_dir

cmd="$time -v $star --runMode genomeGenerate --runThreadN $n_threads --genomeDir $starsim_star_idx_dir --genomeFastaFiles $humanCR3_genome_path --sjdbGTFfile $humanCR3_genes_path > $starsim_star_idx_dir/star_genomeGenerate.time 2>&1"
eval $cmd

# quantification
starsim_star_quant_dir="$starsim_star_dir/star_quant"
mkdir -p $starsim_star_quant_dir

# prepare whitelist
whitelist_path="$ALEVIN_FRY_HOME/plist/10x_v3_permit.txt"

cmd="/usr/bin/time -v $star --genomeDir $starsim_star_idx_dir --soloFeatures Gene --runThreadN $n_threads --readFilesIn $starsim_read2_path $starsim_read1_path --soloCBwhitelist $whitelist_path --soloUMIlen 12 --limitIObufferSize 50000000 50000000 --soloType CB_UMI_Simple --outSAMtype None --outFileNamePrefix ${starsim_star_quant_dir}/ > $starsim_star_quant_dir/starsolo.time 2>&1"

eval $cmd

# gzip quant result to ease the analysis
cmd="gzip ${starsim_star_quant_dir}/Solo.out/Gene/raw/*"
echo $cmd

cmd="gzip ${starsim_star_quant_dir}/Solo.out/Gene/filtered/*"
echo $cmd

#---------------------------------------------------------------------------------------------------------------#
echo "  - Running kallisto-D|bustools"
# Kallisto-D|bustools
## define working dir
starsim_kb_dir="$starsim_root_dir/kallisto_bustools_results"
mkdir -p $starsim_kb_dir

## Build index kallisto|bustools with D-list
### make dir
cd $starsim_kb_dir
starsim_kb_idx_dir="$starsim_kb_dir/kb_index"
mkdir -p $starsim_kb_idx_dir
starsim_kb_single_cell_idx_path="$starsim_kb_idx_dir/kbd_idx_mature_as_ref_genome_as_dlist.idx"

#### kb ref is used to extract the matrure transcripts
cmd="$time -v $kb ref -i $starsim_kb_idx_dir/standard_index.idx --kallisto $kallistod --workflow standard --overwrite -f1 $starsim_kb_idx_dir/f1 -g $starsim_kb_idx_dir/g $humanCR3_genome_path  $humanCR3_genes_path > $starsim_kb_idx_dir/kb_ref.time 2>&1"
eval $cmd

#### for cytoplamic (single-cell) dataset, the index is constructed using
#### mature transcripts as the reference, and the whole genome as the D-list
cmd="$time -v $kallistod index -t $n_threads -i $starsim_kb_single_cell_idx_path -d $humanCR3_genome_path $starsim_kb_idx_dir/f1 > $starsim_kb_idx_dir/kbd_idx_mature_as_ref_genome_as_dlist.time 2>&1"
eval $cmd

#### quantification
starsim_kb_quant_dir="$starsim_kb_dir/kb_quant"
mkdir -p $starsim_kb_quant_dir

cmd="$time -v $kb count --overwrite --kallisto $kallistod --bustools $bustools --workflow standard -i $starsim_kb_single_cell_idx_path -g $starsim_kb_idx_dir/g -t $n_threads -x 10XV3 -o $starsim_kb_quant_dir $starsim_read1_path $starsim_read2_path > $starsim_kb_quant_dir/kb_count.time 2>&1"
eval $cmd