#!/usr/bin/env bash

#  This script serves as an entry point of the experiments performed in the manuscript.
## There are totally three experiemnts:
## 1. Read splicing status classification experiment
## 2. Analysis of a mouse single-nucleus RNA-seq dataset
## (optional) Analysis of STARsolo simulation

script_dir="$(dirname "$0")"
echo $script_dir

#################################################################################################################
# classification experiment
#################################################################################################################
# There are three parts in this experiment
# 1. read simulation
# 2. read mapping
# 3. evaluation
echo "================================================================="
echo "Start running the classification experiment"
sh $script_dir/run_classification_experiment.sh $script_dir



#################################################################################################################
# Real mouse neuron nuclei dataset
#################################################################################################################
# As the dataset referred in the preprint and the GitHub repo are not the same  
# In the preprint: https://www.10xgenomics.com/resources/datasets/5k-adult-mouse-brain-nuclei-isolated-with-chromium-nuclei-isolation-kit-3-1-standard
# in the Github Repo: https://www.10xgenomics.com/resources/datasets/5-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-nuclei-3-1-standard-6-0-0
# here we qantify both of them  

echo "================================================================="
echo "Start running the mouse nuclei experiment"
sh $script_dir/run_mouse_nuclei.sh $script_dir



#################################################################################################################
# STARsolo simulation
#################################################################################################################
# In this experiment, we run the MultiGeneYes simulation from STARsolo preprint. 
# The quantification result will be generated using all tested methods 
# Then, by comparing with the ground truth, we compute the cell-level Spearman correlation
# of each tested method to the ground truth.

# echo "================================================================="
# echo "Start running the STARsolo simulation experiment"
# sh $script_dir/run_starsolo_simulation.sh $script_dir 
