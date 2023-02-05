# vespa.tutorial - Virtual Enrichment-based Signaling Protein-activity Analysis

## Introduction
This tutorial provides instructions on how to use the the VESPA algorithm, including (1) data import, (2) network reconstruction & (3) activity inference. 

This tutorial processes the [``CPTAC S056 LUAD``](https://pubmed.ncbi.nlm.nih.gov/32649874/) dataset. The original data was obtained from [PDC](https://pdc.cancer.gov) and is licensed under [Creative Commons CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) licensing terms. Please contact the original authors and PDC for any queries regarding the dataset.

## Installation
This tutorial requires installation of the [``vespa``](https://github.com/califano-lab/vespa) and [``vespa.db``](https://github.com/califano-lab/vespa.db) R-packages. Please check the respective repositories for instructions. Further, [``vespa.net``](https://github.com/califano-lab/vespa.net) and dependencies need to be installed.

### 1. Data import
In this first step, original CPTAC data is imported in CCT file and converted to the VESPA file format:

```
cd 01_import
Rscript import.R
cd ..
```

This step generates two RDS files, ``CPTAC_S046_LUAD_phospho.rds`` and ``CPTAC_S046_LUAD_proteo.rds``.

### 2. Network reconstruction
In the next step, the signaling network is being reconstructed. For this purpose, we copy the files from the import step, as well as the FASTA library to the vespa.net folder:

``
cd 02_network
git clone https://github.com/califano-lab/vespa.net.git
cp ../01_import/CPTAC_S046_LUAD_phospho.rds vespa.net/
cp ../01_import/CPTAC_S046_LUAD_proteo.rds vespa.net/
cp ../01_import/CPTAC_S046_LUAD_phospho.rds vespa.net/reference.rds
cp ../01_import/library.fasta vespa.net/
cd ../
``

Next, we execute the Snakemake workflow (with 8 CPUs). Here we will use Singularity in a HPC environment. This step might need to be adapted to your local HPC or cluster infrastructure. It can also be run on a Linux Workstation with the packages installed locally.

``
cd 02_network/vespa.net
sbatch --qos=1day --time=1-00:00:00 --mem-per-cpu=8192 snakemake --snakefile Snakefile --use-singularity -j 64 --restart-times 2 --cluster-config envs/res.json --cluster "sbatch --ntasks {cluster.nCPUs} --mem-per-cpu {cluster.memory} --qos=6hours --time=6:00:0"
cp results/*.rds ../
cd ../../
``

The output are the signaling networks in RDS format.

### 3. Protein activity inference
In the final step, we will use the phosphoproteomic data and the signaling networks to infer protein activity for the reconstructed kinases and phosphatases using VESPA:

```
cd 03_activity
Rscript vespa.R
cd ..
```

The output are each three CSV and PDF files representing protein activity on subtrate, activity & integrated levels. These data matrices can be used for further analysis & data integration. While data import and network reconstruction typically are universal methods, protein activity inference can be adapted to particular experimental designs analgously to the VIPER algorithm.
