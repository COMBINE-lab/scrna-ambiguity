#!/usr/bin/env bash
script_dir="$1"

. $script_dir/run_me.config


#################################################################################################################
# classification experiment
#################################################################################################################
# There are three parts in this experiment
# 1. read simulation
# 2. read mapping
# 3. evaluation

#---------------------------------------------------------------------------------------------------------------#
# TODO: upload simulated raeds to Zenodo and move the code to a separate script
# define working dir

classification_root_dir="$root_dir/classification_experiment"
mkdir -p $classification_root_dir

classification_data_dir="$classification_root_dir/data"
mkdir -p $classification_data_dir

# download 10x human 2020-A reference
echo " - Fetching 10x human 2020-A reference"
cmd="wget -qO- https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz | tar xzf - -C $classification_data_dir"
eval $cmd

human2020A_ref_dir="$classification_data_dir/refdata-gex-GRCh38-2020-A"
human2020A_genome_path="$human2020A_ref_dir/fasta/genome.fa"
human2020A_genes_path="$human2020A_ref_dir/genes/genes.gtf"

# make spliceu reference
echo " - Making spliceu reference"
human2020A_spliceu_ref_dir="$human2020A_ref_dir/refdata-gex-GRCh38-2020-A_spliceu_ref"
human2020A_spliceu_ref_prefix="refdata-gex-GRCh38-2020-A_spliceu"
mkdir -p $human2020A_spliceu_ref_dir

cmd="Rscript $script_dir/make_spliceu_txome.R $human2020A_genome_path $human2020A_genes_path $human2020A_spliceu_ref_dir $human2020A_spliceu_ref_prefix"
eval $cmd


human2020A_spliceu_mature_path="$human2020A_spliceu_ref_dir/${human2020A_spliceu_ref_prefix}_mature.fa"
human2020A_spliceu_nascent_path="$human2020A_spliceu_ref_dir/${human2020A_spliceu_ref_prefix}_nascent.fa"
human2020A_spliceu_all_path="$human2020A_spliceu_ref_dir/${human2020A_spliceu_ref_prefix}_all.fa"
human2020A_spliceu_t2g_3col_path="$human2020A_spliceu_ref_dir/${human2020A_spliceu_ref_prefix}_t2g_3col.tsv"

# BBmap simulation
echo " - Running BBmap simulation"

bbmap_sim_dir="$classification_data_dir/bbmap_sim"
mkdir -p $bbmap_sim_dir
randomreads="$bbmap/randomreads.sh"

## specify bbmap sim path
bbmap_sim_fastq_dir="$bbmap_sim_dir/read_fastqs"
mkdir $bbmap_sim_fastq_dir

### read1 path
cytoplasmic_mature_R1_path="$bbmap_sim_fastq_dir/cytoplasmic_mature/cytoplasmic_mature_R1.fastq"
cytoplasmic_nascent_R1_path="$bbmap_sim_fastq_dir/cytoplasmic_nascent/cytoplasmic_nascent_R1.fastq"
nuclear_mature_R1_path="$bbmap_sim_fastq_dir/nuclear_mature/nuclear_mature_R1.fastq"
nuclear_nascent_R1_path="$bbmap_sim_fastq_dir/nuclear_nascent/nuclear_nascent_R1.fastq"

### read2 path
cytoplasmic_mature_R2_path="$bbmap_sim_fastq_dir/cytoplasmic_mature/cytoplasmic_mature_R2.fastq"
cytoplasmic_nascent_R2_path="$bbmap_sim_fastq_dir/cytoplasmic_nascent/cytoplasmic_nascent_R2.fastq"
nuclear_mature_R2_path="$bbmap_sim_fastq_dir/nuclear_mature/nuclear_mature_R2.fastq"
nuclear_nascent_R2_path="$bbmap_sim_fastq_dir/nuclear_nascent/nuclear_nascent_R2.fastq"

echo "   - Simulate reads"
cd $bbmap_sim_dir
## simulate reads
### sc mature
cmd="$randomreads seed=42 ref=$human2020A_spliceu_mature_path paired=f out=$cytoplasmic_mature_R2_path q=36 adderrors=f banns=t minlen=150 maxlen=150 reads=4000000"
eval $cmd

### sc nascent
cmd="$randomreads seed=42 ref=$human2020A_spliceu_nascent_path paired=f out=$cytoplasmic_nascent_R2_path q=36 adderrors=f banns=t minlen=150 maxlen=150 reads=1000000"
eval $cmd

### sn mature
cmd="$randomreads seed=42 ref=$human2020A_spliceu_mature_path paired=f out=$nuclear_mature_R2_path q=36 adderrors=f banns=t minlen=150 maxlen=150 reads=1000000"
eval $cmd

### sn nascent
cmd="$randomreads seed=42 ref=$human2020A_spliceu_nascent_path paired=f out=$nuclear_nascent_R2_path q=36 adderrors=f banns=t minlen=150 maxlen=150 reads=4000000"
eval $cmd

## make reads as 10x format
echo "   - Make reads as 10x format"

cmd="Rscript $script_dir/bbmap_to_10x.R $cytoplasmic_mature_R2_path $cytoplasmic_nascent_R2_path $nuclear_mature_R2_path $nuclear_nascent_R2_path $cytoplasmic_mature_R1_path $cytoplasmic_nascent_R1_path $nuclear_mature_R1_path $nuclear_nascent_R1_path $n_threads"
eval $cmd

# get barcode (CB+UMI) to read name mapping
# awk 'BEGIN {capture=0;rname="nothing";rseq="nothing"}{if ( $1 ~ /^\+/ ) {capture=0; print rseq"\t"rname;} else if (capture) {rseq=substr($0,2);} else if ( $1 ~ /^@/ ) {rname=substr($0,2); capture=1;} }'

#---------------------------------------------------------------------------------------------------------------#
echo "   - Get matching based ground truth assignment"

# get matching based ground truth assignment
matching_based_truth_dir="$bbmap_sim_dir/matching_based_truth"
mkdir -p $matching_based_truth_dir

## build FM index
echo "   - Build FM index"
human2020A_spliceu_fm_index_path="$matching_based_truth_dir/human2020A_spliceu_fm_index"
cmd="$empirical_splice_status/build_index $human2020A_spliceu_all_path $human2020A_spliceu_fm_index_path"
eval $cmd

## find occurences
echo "   - Find occurences"

cmd="$empirical_splice_status/map_queries $human2020A_spliceu_fm_index_path $cytoplasmic_mature_R2_path > $matching_based_truth_dir/cytoplasmic_mature_splice_status_truth.tsv"
eval $cmd

cmd="$empirical_splice_status/map_queries $human2020A_spliceu_fm_index_path $cytoplasmic_nascent_R2_path > $matching_based_truth_dir/cytoplasmic_nascent_splice_status_truth.tsv"
eval $cmd

cmd="$empirical_splice_status/map_queries $human2020A_spliceu_fm_index_path $nuclear_mature_R2_path > $matching_based_truth_dir/nuclear_mature_splice_status_truth.tsv"
eval $cmd

cmd="$empirical_splice_status/map_queries $human2020A_spliceu_fm_index_path $nuclear_nascent_R2_path > $matching_based_truth_dir/nuclear_nascent_splice_status_truth.tsv"
eval $cmd

## convert read name in truth files to barcode
echo "   - Convert read name in truth files to barcode"
cmd="Rscript $script_dir/convert_readname_to_barcode.R $matching_based_truth_dir $bbmap_sim_fastq_dir"
eval $cmd

# 2. Read alignment 

#---------------------------------------------------------------------------------------------------------------#
echo " - Running alevin-fry splici"

## alevin-fry
### splici
classification_af_splici_dir="$classification_root_dir/alevin_fry_splici_results"
mkdir -p $classification_af_splici_dir

classification_af_splici_idx_dir="$classification_af_splici_dir/af_splici_idx"
mkdir -p $classification_af_splici_idx_dir

# build index
cmd="$time -v $simpleaf index -o $classification_af_splici_idx_dir -t $n_threads -f $human2020A_genome_path -g $human2020A_genes_path -r 150 > $classification_af_splici_idx_dir/simpleaf_index.time 2>&1"
eval $cmd

classification_af_splici_map_dir="$classification_af_splici_dir/af_splici_map"
mkdir -p $classification_af_splici_map_dir
### mapping
for sname in "cytoplasmic_mature" "cytoplasmic_nascent" "nuclear_mature" "nuclear_nascent"
do
    o_dir="$classification_af_splici_map_dir/$sname"
    mkdir -p $o_dir

    cmd="$time -v $salmon alevin --chromiumV3 -p $n_threads -i $classification_af_splici_idx_dir/index -l IU --rad -o $o_dir -1 $bbmap_sim_fastq_dir/$sname/${sname}_R1.fastq -2 $bbmap_sim_fastq_dir/$sname/${sname}_R2.fastq --sketch > $o_dir/salmon_alevin.time 2>&1"
    eval $cmd

    # write mapping record to tsv
    cmd="$af view -r ${o_dir}/map.rad > ${o_dir}/map.rad.tsv"
    eval $cmd

    # get the final splice status
    cmd="Rscript $script_dir/empirical_splice_status_af.R ${o_dir}/map.rad.tsv  ${classification_af_splici_map_dir}/${sname}_splice_status_empirical.tsv"
    eval $cmd
done

#---------------------------------------------------------------------------------------------------------------#
echo " - Running alevin-fry spliceu"
### spliceu
# spliceu is here $human2020A_spliceu_all_path

classification_af_spliceu_dir="$classification_root_dir/alevin_fry_spliceu_results"
mkdir -p $classification_af_spliceu_dir

classification_af_spliceu_idx_dir="$classification_af_spliceu_dir/af_spliceu_idx"
mkdir -p $classification_af_spliceu_idx_dir

# build index
cmd="$time -v $simpleaf index --keep-duplicates -o $classification_af_spliceu_idx_dir -t $n_threads --ref-seq $human2020A_spliceu_all_path > $classification_af_spliceu_idx_dir/simpleaf_index.time 2>&1"
eval $cmd

classification_af_spliceu_map_dir="$classification_af_spliceu_dir/af_spliceu_map"
mkdir -p $classification_af_spliceu_map_dir
### mapping
for sname in "cytoplasmic_mature" "cytoplasmic_nascent" "nuclear_mature" "nuclear_nascent"
do
    o_dir="$classification_af_spliceu_map_dir/$sname"
    mkdir -p $o_dir

    cmd="$time -v $salmon alevin --chromiumV3 -p $n_threads -i $classification_af_spliceu_idx_dir/index -l IU --rad -o $o_dir -1 $bbmap_sim_fastq_dir/$sname/${sname}_R1.fastq -2 $bbmap_sim_fastq_dir/$sname/${sname}_R2.fastq --sketch > $o_dir/salmon_alevin.time 2>&1"
    eval $cmd

    # write mapping record to tsv
    cmd="$af view -r ${o_dir}/map.rad > ${o_dir}/map.rad.tsv"
    eval $cmd

    # get the final splice status
    cmd="Rscript $script_dir/empirical_splice_status_af.R ${o_dir}/map.rad.tsv  ${classification_af_spliceu_map_dir}/${sname}_splice_status_empirical.tsv"
    eval $cmd
done

#---------------------------------------------------------------------------------------------------------------#
echo " - Running alevin-fry spliceu piscem"
### spliceu piscem
classification_af_spliceu_piscem_dir="$classification_root_dir/alevin_fry_spliceu_piscem_results"
mkdir -p $classification_af_spliceu_piscem_dir

# build index
classification_af_spliceu_piscem_idx_dir="$classification_af_spliceu_piscem_dir/af_spliceu_piscem_idx"
mkdir -p $classification_af_spliceu_piscem_idx_dir

cmd="$time -v $piscem build --overwrite -s $human2020A_spliceu_all_path -k 31 -m 19 -t $n_threads -o $classification_af_spliceu_piscem_idx_dir/af_spliceu_piscem > $classification_af_spliceu_piscem_idx_dir/piscem_build.time 2>&1"
eval $cmd

### mapping
classification_af_spliceu_piscem_map_dir="$classification_af_spliceu_piscem_dir/piscem_map"
mkdir -p $classification_af_spliceu_piscem_map_dir

for sname in "cytoplasmic_mature" "cytoplasmic_nascent" "nuclear_mature" "nuclear_nascent"
do
    o_dir="$classification_af_spliceu_piscem_map_dir/$sname"
    mkdir -p $o_dir

    cmd="$time -v $piscem map-sc -i $classification_af_spliceu_piscem_idx_dir/af_spliceu_piscem -g chromium_v3 -1 $bbmap_sim_fastq_dir/$sname/${sname}_R1.fastq -2 $bbmap_sim_fastq_dir/$sname/${sname}_R2.fastq -t $n_threads -o $o_dir > $o_dir/piscem_map_sc.time 2>&1"
    eval $cmd

    # write mapping record to tsv
    cmd="$af view -r ${o_dir}/map.rad > ${o_dir}/map.rad.tsv"
    eval $cmd
    
    # get the final splice status
    cmd="Rscript $script_dir/empirical_splice_status_af.R ${o_dir}/map.rad.tsv  ${classification_af_spliceu_piscem_map_dir}/${sname}_splice_status_empirical.tsv"
    eval $cmd
done

#---------------------------------------------------------------------------------------------------------------#
echo " - Running kallisto-D|bustools"
## Kallisto-D|bustools
### define working dir
classification_kb_dir="$classification_root_dir/kallisto_bustools_results"
mkdir -p $classification_kb_dir

### Build index
### kallisto|bustools with D-list
#### make dir
classification_kb_idx_dir="$classification_kb_dir/kb_index"
mkdir -p $classification_kb_idx_dir
classification_kb_single_cell_idx_path="$classification_kb_idx_dir/kbd_idx_mature_as_ref_genome_as_dlist.idx"
classification_kb_single_nucleus_idx_path="$classification_kb_idx_dir/kbd_idx_nascent_as_ref_mature_as_dlist.idx"

##### kb ref is used to extract the matrure transcripts
cmd="$time -v $kb ref -i $classification_kb_idx_dir/standard_index.idx --kallisto $kallistod --workflow standard --overwrite -f1 $classification_kb_idx_dir/f1 -g $classification_kb_idx_dir/g $human2020A_genome_path  $human2020A_genes_path > $classification_kb_idx_dir/kb_ref.time 2>&1"
eval $cmd

##### for cytoplamic (single-cell) dataset, the index is constructed using
##### mature transcripts as the reference, and the whole genome as the D-list
cmd="$time -v $kallistod index -t $n_threads -i $classification_kb_single_cell_idx_path -d $human2020A_genome_path $classification_kb_idx_dir/f1 > $classification_kb_idx_dir/kbd_idx_mature_as_ref_genome_as_dlist.time 2>&1"
eval $cmd

##### for nuclear (single-nucleus) dataset, the index is constructed using
##### nascent transcripts as the reference, and the mature transcripts as the D-list
cmd="$time -v $kallistod index -t $n_threads -i $classification_kb_single_nucleus_idx_path -d $classification_kb_idx_dir/f1 $human2020A_spliceu_nascent_path > $classification_kb_idx_dir/kbd_idx_nascent_as_ref_mature_as_dlist.time 2>&1"
eval $cmd

### mapping
#### make directory
classification_kb_map_dir="$classification_kb_dir/kbd_map"

for sname in "cytoplasmic_mature" "cytoplasmic_nascent" "nuclear_mature" "nuclear_nascent"
do
    for index in 'kbd_idx_mature_as_ref_genome_as_dlist' 'kbd_idx_nascent_as_ref_mature_as_dlist'
    do
        # map
        o_dir="$classification_kb_map_dir/$sname/$index"
        mkdir -p $o_dir

        cmd="$time -v $kallistod bus -x 10xv3 --unstranded -i "$classification_kb_idx_dir/${index}.idx" -o $o_dir -t $n_threads $bbmap_sim_fastq_dir/$sname/${sname}_R1.fastq $bbmap_sim_fastq_dir/$sname/${sname}_R2.fastq > $o_dir/kallisto_bus.time 2>&1"
        eval $cmd

        # write mapping record to tsv
        cmd="$bustools text $o_dir/output.bus --output $o_dir/output.bus.tsv"
        eval $cmd
    done

    # after mapping a sample to both indices, get the final splice status
    cmd="Rscript $script_dir/empirical_splice_status_kb.R $classification_kb_map_dir/$sname/kbd_idx_mature_as_ref_genome_as_dlist/output.bus.tsv $classification_kb_map_dir/$sname/kbd_idx_nascent_as_ref_mature_as_dlist/output.bus.tsv $classification_kb_map_dir/${sname}_splice_status_empirical.tsv"
    eval $cmd
done

