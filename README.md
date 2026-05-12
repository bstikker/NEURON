# 📦 NeuroTide

(add your badges here)

> *Your documentation is a direct reflection of your software, so hold it to the same standards.*


## 🌟 Highlights

- Some functionality made easy!
- This problem handled
- etc.


## ℹ️ Overview

This snake pipeline was developed for the molecular diagnostics department of Amsterdam UMC and can be used for automated analysis of nanopore-sequencing data, used for diagnostics of adult-type diffuse gliomas.
It is tailored for the SLURM-controlled HPC environment 'Helios' of Amsterdam UMC, and containerized for Singularity/docker use. Starting from basecalled and aligned nanopore bam files, the pipeline includes the Sturgeon v1 and Sturgeon v2 CNS tumor classifiers (Vermeulen et al. 2023), MGMT promoter methylation status (Patel et al. 2022), CNV profiling using QDNAseq, and tumor cell fraction estimation infered from the copy number data using ACE (ref), and compiles it nicely into 1 PDF report.
In addition, this pipeline also offers a reconstructed timeline analysis based on the ONT's sequencing summary file. You can predefine the time bins. Each bin gives you Sturgeon V1 and V2 outcomes as well as a bin-corresponding CNV profile. This is visualized in line graphs over time, and again compiled into 1 PDF report.


### ✍️ Authors

Bernard Stikker - Department of Pathology - Amsterdam UMC


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

First install snakemake
Install snakemake SLURM plugin
Install singularity
Install from .def file


```bash
pip install my-package
```

And be sure to specify any other minimum requirements like Python versions or operating systems.

*You may be inclined to add development instructions here, don't.*


## 💭 Feedback and Contributing

Add a link to the Discussions tab in your repo and invite users to open issues for bugs/feature requests.

This is also a great place to invite others to contribute in any ways that make sense for your project. Point people to your DEVELOPMENT and/or CONTRIBUTING guides if you have them.
