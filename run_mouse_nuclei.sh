#!/usr/bin/env bash

script_dir="$1"

. $script_dir/run_me.config

#################################################################################################################
# Real mouse neuron nuclei dataset
#################################################################################################################
# As the dataset referred in the preprint and the GitHub repo are not the same  
# In the preprint: https://www.10xgenomics.com/resources/datasets/5k-adult-mouse-brain-nuclei-isolated-with-chromium-nuclei-isolation-kit-3-1-standard
# in the Github Repo: https://www.10xgenomics.com/resources/datasets/5-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-nuclei-3-1-standard-6-0-0
# here we qantify both of them  

# make working dir
realdata_data_dir="$brain_root_dir/data"
mkdir -p $realdata_data_dir
#---------------------------------------------------------------------------------------------------------------#
# prepare data
# fastq
read_length=90 # from the webpage
realdata_adult_brain_read_dir="$realdata_data_dir/adult_brain_read_fastqs"
realdata_E18_brain_read_dir="$realdata_data_dir/E18_brain_read_fastqs"
mkdir -p $realdata_adult_brain_read_dir
mkdir -p $realdata_E18_brain_read_dir

# we need to use find command to get the file names for each tool later
echo "  - Fetching read fastq files"
cmd="wget -qO- https://cf.10xgenomics.com/samples/cell-exp/7.0.0/5k_mouse_brain_CNIK_3pv3/5k_mouse_brain_CNIK_3pv3_fastqs.tar | tar xf - -C $realdata_adult_brain_read_dir"
eval $cmd

cmd="wget -qO- https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_Nuclei_5K_Multiplex/SC3_v3_NextGem_DI_Nuclei_5K_Multiplex_fastqs.tar | tar xf - -C $realdata_E18_brain_read_dir"
eval $cmd

## make read input for each tool
### alevin-fry
### For simpleaf, fastq files need to be comma separated 
realdata_af_adult_brain_read1_path="$(find -L ${realdata_adult_brain_read_dir} -name "*_R1_*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd,)"
realdata_af_adult_brain_read2_path="$(find -L ${realdata_adult_brain_read_dir} -name "*_R2_*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd,)"
realdata_af_E18_brain_read1_path="$(find -L ${realdata_E18_brain_read_dir} -name "*_R1_*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd,)"
realdata_af_E18_brain_read2_path="$(find -L ${realdata_E18_brain_read_dir} -name "*_R2_*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd,)"

### starsolo
realdata_ss_adult_brain_read_path="$realdata_af_adult_brain_read2_path $realdata_af_adult_brain_read1_path"
realdata_ss_E18_brain_read_path="$realdata_af_E18_brain_read2_path $realdata_af_E18_brain_read1_path"

# kallisto|bustools
realdata_kb_adult_brain_read_path="$(find -L $realdata_adult_brain_read_dir -name "*_R*" -type f | sort | awk '{$1=$1;print}' | paste -sd' ')"
realdata_kb_E18_brain_read_path="$(find -L $realdata_E18_brain_read_dir -name "*_R*" -type f | sort | awk '{$1=$1;print}' | paste -sd' ')"

# reference
echo "  - Fetching Mus musculus GRCm39 108 reference"

mouseGRCm39_ref_dir="$realdata_data_dir/Mus_musculus.GRCm39.108"
mkdir -p  $mouseGRCm39_ref_dir

mouseGRCm39_genome_path="$mouseGRCm39_ref_dir/Mus_musculus.GRCm39.dna.primary_assembly.fa"
mouseGRCm39_genes_path="$mouseGRCm39_ref_dir/Mus_musculus.GRCm39.108.gtf"

cmd="wget -qO- https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz | gunzip - > $mouseGRCm39_genome_path"
eval $cmd

cmd="wget -qO- https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz | gunzip - > $mouseGRCm39_genes_path"
eval $cmd

# make spliceu reference
echo "  - Making spliceu reference"
mouseGRCm39_spliceu_ref_dir="$mouseGRCm39_ref_dir/Mus_musculus.GRCm39.108_spliceu_ref"
mouseGRCm39_spliceu_ref_prefix="Mus_musculus.GRCm39.108_spliceu"
mkdir -p $mouseGRCm39_spliceu_ref_dir

cmd="$time -v Rscript $script_dir/make_spliceu_txome.R $mouseGRCm39_genome_path $mouseGRCm39_genes_path $mouseGRCm39_spliceu_ref_dir $mouseGRCm39_spliceu_ref_prefix > $mouseGRCm39_spliceu_ref_dir/make_spliceu_txome.time 2>&1"
eval $cmd

mouseGRCm39_spliceu_mature_path="$mouseGRCm39_spliceu_ref_dir/${mouseGRCm39_spliceu_ref_prefix}_mature.fa"
mouseGRCm39_spliceu_nascent_path="$mouseGRCm39_spliceu_ref_dir/${mouseGRCm39_spliceu_ref_prefix}_nascent.fa"
mouseGRCm39_spliceu_all_path="$mouseGRCm39_spliceu_ref_dir/${mouseGRCm39_spliceu_ref_prefix}_all.fa"
mouseGRCm39_spliceu_t2g_3col_path="$mouseGRCm39_spliceu_ref_dir/${mouseGRCm39_spliceu_ref_prefix}_t2g_3col.tsv"
mouseGRCm39_spliceu_g2g_path="$mouseGRCm39_spliceu_ref_dir/${mouseGRCm39_spliceu_ref_prefix}_g2g.tsv"

# 2. Quantification 

#---------------------------------------------------------------------------------------------------------------#
echo "  - Running alevin-fry splici"

## alevin-fry splici
realdata_af_splici_dir="$brain_root_dir/alevin_fry_splici_results"
mkdir -p $realdata_af_splici_dir

### splici
realdata_af_splici_idx_dir="$realdata_af_splici_dir/af_splici_idx"
mkdir -p $realdata_af_splici_idx_dir

# build index
cmd="$time -v $simpleaf index -o $realdata_af_splici_idx_dir -t $n_threads -f $mouseGRCm39_genome_path -g $mouseGRCm39_genes_path -r $read_length > $realdata_af_splici_idx_dir/simpleaf_index.time 2>&1"
eval $cmd

# quant
### adult mouse brain nuclei dataset
realdata_adult_brain_af_splici_quant_dir="$realdata_af_splici_dir/af_splici_quant/adult_brain"
mkdir -p $realdata_adult_brain_af_splici_quant_dir

cmd="$time -v $simpleaf quant -c 10xv3 -o $realdata_adult_brain_af_splici_quant_dir -t $n_threads -i $realdata_af_splici_idx_dir/index  -u -r cr-like -m $realdata_af_splici_idx_dir/index/t2g_3col.tsv -1 $realdata_af_adult_brain_read1_path -2 $realdata_af_adult_brain_read2_path > $realdata_adult_brain_af_splici_quant_dir/simpleaf_quant.time 2>&1"
eval $cmd

### combined Cortex, Hippocampus and Subventricular Zone Nuclei
realdata_E18_brain_af_splici_quant_dir="$realdata_af_splici_dir/af_splici_quant/E18_brain"
mkdir -p $realdata_E18_brain_af_splici_quant_dir

cmd="$time -v $simpleaf quant -c 10xv3 -o $realdata_E18_brain_af_splici_quant_dir -t $n_threads -i $realdata_af_splici_idx_dir/index  -u -r cr-like -m $realdata_af_splici_idx_dir/index/t2g_3col.tsv -1 $realdata_af_E18_brain_read1_path -2 $realdata_af_E18_brain_read2_path > $realdata_E18_brain_af_splici_quant_dir/simpleaf_quant.time 2>&1"
eval $cmd

#---------------------------------------------------------------------------------------------------------------#
echo "  - Running alevin-fry spliceu"
### spliceu
realdata_af_spliceu_dir="$brain_root_dir/alevin_fry_spliceu_results"

mkdir -p $realdata_af_spliceu_dir

realdata_af_spliceu_idx_dir="$realdata_af_spliceu_dir/af_spliceu_idx"
mkdir -p $realdata_af_spliceu_idx_dir

# build index
cmd="$time -v $simpleaf index --keep-duplicates -o $realdata_af_spliceu_idx_dir -t $n_threads --ref-seq $mouseGRCm39_spliceu_all_path > $realdata_af_spliceu_idx_dir/simpleaf_index.time 2>&1"
eval $cmd

# quant
## adult brain
realdata_adult_brain_af_spliceu_quant_dir="$realdata_af_spliceu_dir/af_spliceu_quant/adult_brain"
mkdir -p $realdata_adult_brain_af_spliceu_quant_dir
cmd="$time -v $simpleaf quant -c 10xv3 -o $realdata_adult_brain_af_spliceu_quant_dir -t $n_threads -i $realdata_af_spliceu_idx_dir/index  -u -r cr-like -m $mouseGRCm39_spliceu_t2g_3col_path -1 $realdata_af_adult_brain_read1_path -2 $realdata_af_adult_brain_read2_path > $realdata_adult_brain_af_spliceu_quant_dir/simpleaf_quant.time 2>&1"
eval $cmd

### combined Cortex, Hippocampus and Subventricular Zone Nuclei
realdata_E18_brain_af_spliceu_quant_dir="$realdata_af_spliceu_dir/af_spliceu_quant/E18_brain"
mkdir -p $realdata_E18_brain_af_spliceu_quant_dir
cmd="$time -v $simpleaf quant -c 10xv3 -o $realdata_E18_brain_af_spliceu_quant_dir -t $n_threads -i $realdata_af_spliceu_idx_dir/index  -u -r cr-like -m $mouseGRCm39_spliceu_t2g_3col_path -1 $realdata_af_E18_brain_read1_path -2 $realdata_af_E18_brain_read2_path > $realdata_E18_brain_af_spliceu_quant_dir/simpleaf_quant.time 2>&1"
eval $cmd


#---------------------------------------------------------------------------------------------------------------#
echo "  - Running alevin-fry spliceu piscem"
# spliceu piscem
realdata_af_spliceu_piscem_dir="$brain_root_dir/alevin_fry_spliceu_piscem_results"
mkdir -p $realdata_af_spliceu_piscem_dir

# build index
realdata_af_spliceu_piscem_idx_dir="$realdata_af_spliceu_piscem_dir/af_spliceu_piscem_idx"
mkdir -p $realdata_af_spliceu_piscem_idx_dir

cmd="$time -v $piscem build -s $mouseGRCm39_spliceu_all_path -k 31 -m 19 -t $n_threads -o $realdata_af_spliceu_piscem_idx_dir/af_spliceu_piscem > $realdata_af_spliceu_piscem_idx_dir/piscem_build.time 2>&1"
eval $cmd

# quant
## adult brain
realdata_adult_brain_af_spliceu_piscem_quant_dir="$realdata_af_spliceu_piscem_dir/af_spliceu_piscem_quant/adult_brain"
mkdir -p $realdata_adult_brain_af_spliceu_piscem_quant_dir

# mapping
realdata_adult_brain_af_spliceu_piscem_map_dir="$realdata_adult_brain_af_spliceu_piscem_quant_dir/piscem_map"
mkdir -p $realdata_adult_brain_af_spliceu_piscem_map_dir

cmd="$time -v $piscem map-sc -i $realdata_af_spliceu_piscem_idx_dir/af_spliceu_piscem -g chromium_v3 -1 $realdata_af_adult_brain_read1_path -2 $realdata_af_adult_brain_read2_path -t $n_threads -o $realdata_adult_brain_af_spliceu_piscem_map_dir > $realdata_adult_brain_af_spliceu_piscem_map_dir/piscem_map_sc.time 2>&1"
eval $cmd

# quantification
cmd="$time -v $simpleaf quant -c 10xv3 -o $realdata_adult_brain_af_spliceu_piscem_quant_dir -t $n_threads --map-dir $realdata_adult_brain_af_spliceu_piscem_map_dir -u -r cr-like -m $mouseGRCm39_spliceu_t2g_3col_path > $realdata_adult_brain_af_spliceu_piscem_quant_dir/simpleaf_quant.time 2>&1"
eval $cmd

### combined Cortex, Hippocampus and Subventricular Zone Nuclei
realdata_E18_brain_af_spliceu_piscem_quant_dir="$realdata_af_spliceu_piscem_dir/af_spliceu_piscem_quant/E18_brain"
mkdir -p $realdata_E18_brain_af_spliceu_piscem_quant_dir

# mapping
realdata_E18_brain_af_spliceu_piscem_map_dir="$realdata_E18_brain_af_spliceu_piscem_quant_dir/piscem_map"
mkdir -p $realdata_E18_brain_af_spliceu_piscem_map_dir

cmd="$time -v $piscem map-sc -i $realdata_af_spliceu_piscem_idx_dir/af_spliceu_piscem -g chromium_v3 -1 $realdata_af_E18_brain_read1_path -2 $realdata_af_E18_brain_read2_path -t $n_threads -o $realdata_E18_brain_af_spliceu_piscem_map_dir > $realdata_E18_brain_af_spliceu_piscem_map_dir/piscem_map_sc.time 2>&1"
eval $cmd

# quantification
cmd="$time -v $simpleaf quant -c 10xv3 -o $realdata_E18_brain_af_spliceu_piscem_quant_dir -t $n_threads --map-dir $realdata_E18_brain_af_spliceu_piscem_map_dir -u -r cr-like -m $mouseGRCm39_spliceu_t2g_3col_path > $realdata_E18_brain_af_spliceu_piscem_quant_dir/simpleaf_quant.time 2>&1"
eval $cmd

#---------------------------------------------------------------------------------------------------------------#
echo "  - Running STARsolo"
# STARsolo
realdata_star_dir="$brain_root_dir/starsolo_results"
mkdir -p $realdata_star_dir

# build index
cd $realdata_star_dir
realdata_star_idx_dir="$realdata_star_dir/star_index"
mkdir -p $realdata_star_idx_dir

cmd="$time -v $star --runMode genomeGenerate --runThreadN $n_threads --genomeDir $realdata_star_idx_dir --genomeFastaFiles $mouseGRCm39_genome_path --sjdbGTFfile $mouseGRCm39_genes_path > $realdata_star_idx_dir/star_genomeGenerate.time 2>&1"
eval $cmd

# quantification
# prepare input
whitelist_path="$ALEVIN_FRY_HOME/plist/10x_v3_permit.txt"

## adult brain
realdata_adult_brain_star_quant_dir="$realdata_star_dir/star_quant/adult_brain"
mkdir -p $realdata_adult_brain_star_quant_dir

cmd="/usr/bin/time -v $star --genomeDir $realdata_star_idx_dir --readFilesCommand zcat --soloFeatures GeneFull --runThreadN $n_threads --readFilesIn $realdata_ss_adult_brain_read_path --soloCBwhitelist $whitelist_path --soloUMIlen 12 --limitIObufferSize 50000000 50000000 --soloType CB_UMI_Simple --outSAMtype None --outFileNamePrefix ${realdata_adult_brain_star_quant_dir}/ > $realdata_adult_brain_star_quant_dir/star_quant.time 2>&1"
eval $cmd

# gzip quant result to ease the analysis
cmd="gzip ${realdata_adult_brain_star_quant_dir}/Solo.out/GeneFull/raw/*"
echo $cmd

cmd="gzip ${realdata_adult_brain_star_quant_dir}/Solo.out/GeneFull/filtered/*"
echo $cmd

### combined Cortex, Hippocampus and Subventricular Zone Nuclei
realdata_E18_brain_star_quant_dir="$realdata_star_dir/star_quant/E18_brain"
mkdir -p $realdata_E18_brain_star_quant_dir

cmd="/usr/bin/time -v $star --genomeDir $realdata_star_idx_dir --readFilesCommand zcat --soloFeatures GeneFull --runThreadN $n_threads --readFilesIn $realdata_ss_E18_brain_read_path --soloCBwhitelist $whitelist_path --soloUMIlen 12 --limitIObufferSize 50000000 50000000 --soloType CB_UMI_Simple --outSAMtype None --outFileNamePrefix ${realdata_E18_brain_star_quant_dir}/ > $realdata_E18_brain_star_quant_dir/star_quant.time 2>&1"
eval $cmd

# gzip quant result to ease the analysis
cmd="gzip ${realdata_E18_brain_star_quant_dir}/Solo.out/GeneFull/raw/*"
echo $cmd

cmd="gzip ${realdata_E18_brain_star_quant_dir}/Solo.out/GeneFull/filtered/*"
echo $cmd

#---------------------------------------------------------------------------------------------------------------#
echo "  - Running kallisto-D|bustools using spliceu nascent transcripts"
# Kallisto-D|bustools using spliceu nascent transcripts
# We will use two versions of nascent transcripts to run kb.
# 1. The nascent transcripts from alevin-fry spliceu reference
# 2. The nascent transcripts generated by the generate_cDNA+introns.py 
#    scrtipt in the HSHMP_2022 GitHub repository

# define working dir
realdata_kb_dir="$brain_root_dir/kallisto_bustools_using_spliceu_nasecent_transcripts_results"
mkdir -p $realdata_kb_dir

# Build index
# kallisto|bustools with D-list
## make dir
cd $realdata_kb_dir
realdata_kb_idx_dir="$realdata_kb_dir/kb_index"
mkdir -p $realdata_kb_idx_dir
realdata_kb_single_nucleus_idx_path="$realdata_kb_idx_dir/kbd_idx_nascent_as_ref_mature_as_dlist.idx"

### kb ref is used to extract the matrure transcripts
cmd="$time -v $kb ref -i $realdata_kb_idx_dir/standard_index.idx --kallisto $kallistod --bustools $bustools --workflow standard --overwrite -f1 $realdata_kb_idx_dir/f1 -g $realdata_kb_idx_dir/g $mouseGRCm39_genome_path  $mouseGRCm39_genes_path > $realdata_kb_idx_dir/kb_ref.time 2>&1"
eval $cmd

### for nuclear (single-nucleus) dataset, the index is constructed using
### nascent transcripts as the reference, and the mature transcripts as the D-list
cmd="$time -v $kallistod index -t $n_threads -i $realdata_kb_single_nucleus_idx_path -d $realdata_kb_idx_dir/f1 $mouseGRCm39_spliceu_nascent_path > $realdata_kb_idx_dir/kbd_idx_nascent_as_ref_mature_as_dlist.time 2>&1"
eval $cmd

### quantification
#### adult brain
realdata_adult_brain_kb_quant_dir="$realdata_kb_dir/kb_quant/adult_brain"
mkdir -p $realdata_adult_brain_kb_quant_dir

cmd="$time -v $kb count --overwrite --kallisto $kallistod --bustools $bustools --workflow standard -i $realdata_kb_single_nucleus_idx_path -g $mouseGRCm39_spliceu_g2g_path -t $n_threads -x 10XV3 -o $realdata_adult_brain_kb_quant_dir $realdata_kb_adult_brain_read_path > $realdata_adult_brain_kb_quant_dir/kb_count.time 2>&1"
eval $cmd

### combined Cortex, Hippocampus and Subventricular Zone Nuclei
realdata_E18_brain_kb_quant_dir="$realdata_kb_dir/kb_quant/E18_brain"
mkdir -p $realdata_E18_brain_kb_quant_dir

cmd="$time -v $kb count --overwrite --kallisto $kallistod --bustools $bustools --workflow standard -i $realdata_kb_single_nucleus_idx_path -g $mouseGRCm39_spliceu_g2g_path -t $n_threads -x 10XV3 -o $realdata_E18_brain_kb_quant_dir $realdata_kb_E18_brain_read_path > $realdata_E18_brain_kb_quant_dir/kb_count.time 2>&1"
eval $cmd

#---------------------------------------------------------------------------------------------------------------#
echo "  - Running kallisto-D|bustools using HSHMP nascent transcripts"
## Kallisto-D|bustools using HSHMP nascent transcripts
## We will use two versions of nascent transcripts to run kb.
## 1. The nascent transcripts from alevin-fry spliceu reference
## 2. The nascent transcripts generated by the generate_cDNA+introns.py 
##    scrtipt in the HSHMP_2022 GitHub repository

### define working dir
realdata_kb_dir="$brain_root_dir/kallisto_bustools_using_HSHMP_nascent_transcripts_results"
mkdir -p $realdata_kb_dir

### Build index
### kallisto|bustools with D-list
#### make dir
cd $realdata_kb_dir
realdata_kb_idx_dir="$realdata_kb_dir/kb_index"
mkdir -p $realdata_kb_idx_dir
realdata_kb_single_nucleus_idx_path="$realdata_kb_idx_dir/kbd_idx_nascent_as_ref_mature_as_dlist.idx"

##### make nascent transcripts
mouseGRCm39_HSHMP_nascent_path="$realdata_kb_idx_dir/${mouseGRCm39_spliceu_ref_prefix}_nascent.fa"

cmd="$HSHMP_2022/extract_introns/generate_cDNA+introns.py --nascent --gtf $mouseGRCm39_genes_path --fa $mouseGRCm39_genome_path --out $mouseGRCm39_HSHMP_nascent_path"
eval $cmd

##### kb ref is used to extract the matrure transcripts
##### We still need this to get the t2g file
cmd="$time -v $kb ref -i $realdata_kb_idx_dir/standard_index.idx --kallisto $kallistod --bustools $bustools --workflow standard --overwrite -f1 $realdata_kb_idx_dir/f1 -g $realdata_kb_idx_dir/g $mouseGRCm39_genome_path  $mouseGRCm39_genes_path > $realdata_kb_idx_dir/kb_ref.time 2>&1"
eval $cmd

##### for nuclear (single-nucleus) dataset, the index is constructed using
##### nascent transcripts as the reference, and the mature transcripts as the D-list
cmd="$time -v $kallistod index -t $n_threads -i $realdata_kb_single_nucleus_idx_path -d $realdata_kb_idx_dir/f1 $mouseGRCm39_HSHMP_nascent_path > $realdata_kb_idx_dir/kbd_idx_nascent_as_ref_mature_as_dlist.time 2>&1"
eval $cmd

##### quantification
### adult brain
realdata_adult_brain_kb_quant_dir="$realdata_kb_dir/kb_quant/adult_brain"
mkdir -p $realdata_adult_brain_kb_quant_dir

cmd="$time -v $kb count --overwrite --kallisto $kallistod --bustools $bustools --workflow standard -i $realdata_kb_single_nucleus_idx_path -g $mouseGRCm39_spliceu_g2g_path -t $n_threads -x 10XV3 -o $realdata_adult_brain_kb_quant_dir $realdata_kb_adult_brain_read_path > $realdata_adult_brain_kb_quant_dir/kb_count.time 2>&1"
eval $cmd

### combined Cortex, Hippocampus and Subventricular Zone Nuclei
realdata_E18_brain_kb_quant_dir="$realdata_kb_dir/kb_quant/E18_brain"
mkdir -p $realdata_E18_brain_kb_quant_dir

cmd="$time -v $kb count --overwrite --kallisto $kallistod --bustools $bustools --workflow standard -i $realdata_kb_single_nucleus_idx_path -g $mouseGRCm39_spliceu_g2g_path -t $n_threads -x 10XV3 -o $realdata_E18_brain_kb_quant_dir $realdata_kb_E18_brain_read_path > $realdata_E18_brain_kb_quant_dir/kb_count.time 2>&1"
eval $cmd
