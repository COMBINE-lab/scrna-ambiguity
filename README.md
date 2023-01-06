# Repository accompanying "Understanding and evaluating ambiguity in single-cell and single-nucleus RNA-sequencing"

This repository accompanies the paper:


> Understanding and evaluating ambiguity in single-cell and single-nucleus RNA-sequencing \
> Dongze He, Charlotte Soneson, Rob Patro \
> bioRxiv 2023.01.04.522742; doi: https://doi.org/10.1101/2023.01.04.522742


**Note**: This repository uses git submodules, please clone it with

```{bash}
git clone --recurse-submodules https://github.com/COMBINE-lab/scrna-ambiguity
cd scrna-ambiguity

```

To run the whole pipeline, please set all the paths in the `run_me.config`, and run the `run_me.sh` in a terminal by calling `. run_me.sh`. Running the classification experiment and the analysis of the two mouse brain nuclei datasets  will take approximately 450 GB of disk space. If you choose to run the analysis of the STARsolo simulation as well, please make sure there is another 200 GB of disk space left on your disk.

## Conda environment setup

### Create conda environment
We recommend run the pipeline in a conda environment. To reduce the time of environment creation, we recommend to use [mamba](https://mamba.readthedocs.io/en/latest/installation.html). 

```sh
mamba create -n scrna-ambiguity r-essentials=4.1 r-doParallel=1.0.17 bioconductor-genomicfeatures=1.46.1 bioconductor-biostrings=2.62.0 bioconductor-bsgenome=1.62.0 r-ggplot2=3.4.0 star=2.7.10b kb-python=0.27.3 simpleaf=0.8.1 -y
conda activate scrna-ambiguity

```

If you want to use conda instead, simply replace the `mamba` in the above command with `conda`. The `conda_env.yml` file in the GitHub repository can be used to create the exact conda env we used to run the pipeline. **Note:** The included conda environment was created under linux, and, in general, the code in the current repository is unlikely to work directly under other operating systems.

## Install stand-alone packages
Some packages are not available on conda, so they need to be manually installed.

### BBMap

[BBMap](https://github.com/BioInfoTools/BBMap) is used to simulate the sequencing reads for the classification experiment. We used BBMap 39.01 in the analysis, which can be directly downloaded from [this](https://versaweb.dl.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz) link. The path to the decompressed folder of the downloaded file should be specified in the `run_me.config` file. To download the path and print the path to the decompressed folder, please run the following bash command.

```sh
wget https://versaweb.dl.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz && tar -xzvf BBMap_39.01.tar.gz
cd bbmap && pwd
cd ..

```

### empirical_splice_status

[`empirical_splice_status`](https://github.com/COMBINE-lab/empirical_splice_status) is used to obtain the matching-based true splicing status of the simulated reads in the classification experiment. `empirical_splice_status` is a rust program. To install it, please make sure rust is [installed](https://www.rust-lang.org/tools/install).


The source code of `empirical_splice_status` will be included in this repository if you follow the above `git clone  --recurse-submodules` command to recursively clone this repository. Then you can run the following command without worrying about the `empirical_splice_status` variable in the `run_me.config` file.

```bash
cargo build --release --manifest-path=empirical_splice_status/Cargo.toml

```

If the repository is not cloned recursively, just clone `empirical_splice_status` before build it. If you build it to somewhere else, please use the path to the `release` folder printed from running the following bash command as the`empirical_splice_status` variable in the `run_me.config`

```sh
rm -r empirical_splice_status
git clone https://github.com/COMBINE-lab/empirical_splice_status.git 
cargo build --release --manifest-path=empirical_splice_status/Cargo.toml
cd empirical_splice_status/target/release
pwd
cd ../../..

```

### piscem
[`piscem`](https://github.com/COMBINE-lab/piscem) is the new latest mapping tools in the alevin-fry ecosystem. If offers more concise and memory frugal index than any other methods tested in this work. To install it, you might need to set your `CXX` and `CC` path. Please use the path to the `piscem` executable printed from running the following bash command as the `piscem` variable in the `run_me.config`.

```sh
git clone https://github.com/COMBINE-lab/piscem.git
cd piscem
cargo build --release
cd target/release && find ${PWD} -name "piscem"
cd ../../..

```

### HSHMP_2022
Please use the printed path from the following bash command as the `HSHMP_2022` variable. In this work, we followed the code in commit [cfd6958be6ab65fc340bf1df5f1b0e77f19ff967](https://github.com/pachterlab/HSHMP_2022/tree/cfd6958be6ab65fc340bf1df5f1b0e77f19ff967).

```sh
git clone https://github.com/pachterlab/HSHMP_2022.git
cd HSHMP_2022
git checkout cfd6958be6ab65fc340bf1df5f1b0e77f19ff967
pwd
cd ..

```

### kallisto-D
Please follow the instructions in the [GitHub repository](https://github.com/pachterlab/kallisto-D/tree/26daa940c25ca1787de8aab8d9c8dda1802006ea) to install it, and use the path to the `kallisto` executable as the `kallistod` variable in the `run_me.config` file. In this work, we used commit [`26daa940c25ca1787de8aab8d9c8dda1802006ea`](https://github.com/pachterlab/kallisto-D/tree/26daa940c25ca1787de8aab8d9c8dda1802006ea).

### bustools
Please follow the instructions in the [GitHub repository](https://github.com/BUStools/bustools/tree/89a8a5ddee4ec1c3bb95dfb42e13f2d2ba86e9a1) to install it, and use the path to the `bustools` executable as the `bustools` variable in the `run_me.config`. In this work, we used the commit [`89a8a5ddee4ec1c3bb95dfb42e13f2d2ba86e9a1`](https://github.com/BUStools/bustools/tree/89a8a5ddee4ec1c3bb95dfb42e13f2d2ba86e9a1).
## Specify the root working directory

The `root_dir` variable in the `run_me.config` file serves as the root working directory of the whole pipeline. After specifying a path, please make sure that there are at least 450GB free space in the disk.

## Run the pipeline

After you have set all the paths in the `run_me.config` file and made sure that there are enough free disk space, the pipeline can be run by calling the `run_me.sh` file. Notice that for all variants of alevin-fry, mapping was performed using pseudoalignment with structural constraints corresponding to the `--sketch` flag in `salmon alevin` (likewise, this is the mapping algorithm implemented in [`piscem`](https://github.com/COMBINE-lab/piscem) so far).

```sh
chmod +x *.sh
chmod +x *.R
. run_me.sh

```

## Generate the plots
The figures used in the manuscript can be generated by following the python notebooks in the `notebooks` directory in this repository. You might need to install the python packages used by the notebooks.
- The `classification_experiment.ipynb` is used to generate Figure 2 and 3.
- The `sn_mouse_analysis_adult.ipynb` and `sn_mouse_analysis_E18.ipynb` notebooks are used to generate Figure 4, 5, 6, 7, and 8. Notice that the `sn_mouse_analysis_adult.ipynb` notebook has to be run before `sn_mouse_analysis_E18.ipynb`, as it outputs a file the other one runs.

Notice that the filtered cellular barcode list of the [E18 mouse brain nuclei dataset](https://www.10xgenomics.com/resources/datasets/5-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-nuclei-3-1-standard-6-0-0) can be downloaded from [this link](https://umd.box.com/shared/static/teohgo005zmoq3ha2zxrzviffduvre7w.tsv). As the filtered count matrix is not provided in the dataset's webpage, the filtered barcodes are obtained by loading the `Loupe Browser visualization file (CLOUPE)` file from the dataset's webpage in the Loupe browser, and exporting the clustering. In the final barcode list file, only the barcodes are included, not the clusters. 
For [the adult mouse nuclei dataset](https://www.10xgenomics.com/resources/datasets/5k-adult-mouse-brain-nuclei-isolated-with-chromium-nuclei-isolation-kit-3-1-standard), the filtered barcode list TSV file can be downloaded from [this link](https://umd.box.com/shared/static/p3v83a6oxgmfw1lcwaufkfe3dvr40qjg.tsv). As the filtered barcode list is included in the filtered count matrix file, `Feature / cell matrix (filtered)`, the barcode included in the filtered count matrix is used as the filtered barcode list.

### (OPTIONAL) Run the Analysis of the STARsolo simulation

To run the analysis of the STARsolo simulation, please comment out the line #47-49 in the `run_me.sh`. **As this experiment is not mentioned in the manuscript, this analysis is optional.** 












