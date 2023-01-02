# scrna-ambiguity

**Note**: This repository uses git submodules, please clone it with

```{bash}
git clone --recurse-submodules https://github.com/COMBINE-lab/scrna-ambiguity
```

To run the whole pipeline, please set all the paths in the `run_me.config`, and run the `run_me.sh` in a terminal by calling `. run_me.sh`. Running the classification experiment and the analysis of the two mouse brain nuclei datasets  will take approximately 450 GB of disk space. If you choose to run the analysis of the STARsolo simulation as well, please make sure there is another 200 GB of disk space left on your disk.

## Conda environment setup

### Create conda environment
We recommend run the pipeline in a conda environment. To reduce the time of environment creation, we recommend to use [mamba](https://mamba.readthedocs.io/en/latest/installation.html). 

```sh
mamba create -n scrna-ambiguity r-essentials=4.1 r-doParallel=1.0.17 bioconductor-genomicfeatures=1.46.1 bioconductor-biostrings=2.62.0 bioconductor-bsgenome bsgenome=1.62.0 r-ggplot2=3.4.0 star=2.7.10b kb-python=0.27.3 simpleaf=0.7.0 -y && conda activate scrna-ambiguity
```

If you want to use conda instead, simply replace the `mamba` in the above command with `conda`. The `conda_env.yml` file in the GitHub repository can be used to create the exact conda env we used to run the pipeline. **Note:** The included conda environment was created under linux, and, in general, the code in the current repository is unlikely to work directly under other operating systems.

## Install stand-alone packages
Some packages are not available on conda, so they need to be manually installed.

### BBMap

[BBMap](https://github.com/BioInfoTools/BBMap) is used to simulate the sequencing reads for the classification experiment. We used BBMap 39.01 in the analysis, which can be directly downloaded from [this](https://versaweb.dl.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz) link. The path to the decompressed folder of the downloaded file should be specified in the `run_me.config` file. To download the path and print the path to the decompressed folder, please run the following bash command.

```sh
wget https://versaweb.dl.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz && tar -xzvf BBMap_39.01.tar.gz

cd bbmap && pwd

```

### empirical_splice_status

[`empirical_splice_status`](https://github.com/COMBINE-lab/empirical_splice_status) is used to obtain the matching-based true splicing status of the simulated reads in the classification experiment. `empirical_splice_status` is a rust program. To install it, please make sure rust is [installed](https://www.rust-lang.org/tools/install), and then run the following bash command. Please use the path to the `release` folder printed from running the following bash command as the`empirical_splice_status` variable in the `run_me.config`

```sh
git clone https://github.com/COMBINE-lab/empirical_splice_status.git && cd empirical_splice_status
cargo build --release
cd target/release && pwd

```

### piscem
[`piscem`](https://github.com/COMBINE-lab/piscem) is the new latest mapping tools in the alevin-fry ecosystem. If offers more concise and memory frugal index than any other methods tested in this work. To install it, you might need to set your `CXX` and `CC` path. Please use the path to the `piscem` executable printed from running the following bash command as the `piscem` variable in the `run_me.config`.

```sh
git clone https://github.com/COMBINE-lab/piscem.git && cd piscem
cargo build --release

cd target/release && find ${PWD} -name "piscem"
```

### HSHMP_2022
Please use the printed path from the following bash command as the `HSHMP_2022` variable.

```sh
git clone https://github.com/pachterlab/HSHMP_2022.git && cd HSHMP_2022 && git checkout 03456b623f5c2bb12212b4745e3523cbba57b44c
pwd

```

### kallisto-D
Please follow the instructions in the [GitHub repository](https://github.com/pachterlab/kallisto-D) to install it, and use the path to the `kallisto` executable as the `kallistod` variable in the `run_me.config` file. 


## Specify the root working directory

The `root_dir` variable in the `run_me.config` file serves as the root working directory of the whole pipeline. After specifying a path, please make sure that there are at least 450GB free space in the disk.

## Run the pipeline

After you have set all the paths in the `run_me.config` file and made sure that there are enough free disk space, the pipeline can be run by calling the `run_me.sh` file. 

```sh
chmod +x *.sh
chmod +x *.R
. run_me.sh
```

### (OPTIONAL) Run the Analysis of the STARsolo simulation

To run the analysis of the STARsolo simulation, please comment out the line #65 in the `run_me.sh`. **As this experiment is not mentioned in the manuscript, this analysis is optional.** 










