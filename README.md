# 📦 NEURON

## 🌟 Highlights

- Nanopore methylation-based brain tumor diagnosis (Sturgeon)
- Copy Number Profiles (QDNAseq)
- Estimated tumor purity (ACE)
- MGMT promoter methylation status (Rapid-CNS2)


## ℹ️ Overview

NEURON (NEURO-Nanopore) was developed for the molecular diagnostics department of Amsterdam UMC and can be used for automated analysis of nanopore-sequencing data, used for diagnostics of brain tumors
It is tailored for the SLURM-controlled HPC environment 'Helios' of Amsterdam UMC, and containerized for Singularity/docker use. Starting from basecalled and aligned nanopore bam files, the pipeline includes the Sturgeon v1 and Sturgeon v2 CNS tumor classifiers (Vermeulen et al. 2023), MGMT promoter methylation status (Patel et al. 2022), CNV profiling using QDNAseq, and tumor cell fraction estimation infered from the copy number data using ACE (Poell et al. 2019), and compiles it into 1 PDF report.

In addition, this pipeline also offers a reconstructed timeline analysis based on the ONT's sequencing summary file. You can predefine the time bins. Each bin gives you Sturgeon V1 and V2 outcomes as well as a bin-corresponding CNV profile. This is visualized in line graphs over time, and again compiled into 1 PDF report.

<img width="3466" height="2100" alt="NeuroTide pipeline" src="https://github.com/user-attachments/assets/227527e3-0b24-46fc-b521-eb6dfc5ac934" />

### ✍️ Authors

Bernard Stikker, PhD - Department of Pathology - Amsterdam UMC


## 🚀 Usage

First, we define our sample names. The pipeline will look in the input/ folder for bam files and it's relative location should be defined in the config/samples.yaml file.
ADD EXAMPLE. 

You can configure your SLURM settings in the profile/config.yaml file. In this way, you can run the pipeline as follows:

```bash
snakemake --profile profile 
```
In addition, the pipeline also offers a timecourse analysis. In your config/samples.yaml add the location of ONT's sequencing summary file as well as the name you gave your samples prior to the sequencing run.
This feature can be run as follows:

```bash
snakemake --profile profile rule_timecourse
```

## ⬇️ Installation

#If you haven't installed snakemake yet, this pipeline has been validated using version 9.8.1.
```bash
conda install -c conda-forge -c bioconda snakemake=9.8.1
```
#As Helios uses SLURM as task manager, install the SLURM plugin
```bash
pip install snakemake-executor-plugin-slurm
```
#Install singularity
Helios has singularity installed as module (1.4.1-1.el9)

#Install container from .def file
```bash
apptainer build --fakeroot cns_AUMC_pipeline_v2.2.sif cns_AUMC_pipeline_v2.2.def
```
#Downloading the classifiers (software is installed through the .def file, but models need to be downloaded manually)
Sturgeon V1: https://github.com/UMCUGenetics/sturgeon
Sturgeon V2: Yet to be published

Move the model zip files to references/models/

## Important notes
At the time of writing (26-06-2026) Sturgeon V2 is yet to be published. So V2 functionality is not available in this repo. One published, the singularity container definitions file will be updated to include automated installation. The pipeline currently uses a pre-built container that includes a local copy of the V2 classifier.

## 💭 Feedback and Contributing

Add a link to the Discussions tab in your repo and invite users to open issues for bugs/feature requests.
Please do not hesitate to report bugs or to discuss issues!: https://github.com/bstikker/NEURON/issues
