# Portable pipeline notes

Place these large reference files locally after cloning the repository:

- reference/models/general.zip
- reference/genome/chm13v2.0.fa.gz
- singularity/cns_AUMC_pipeline_v2.sif (build once from singularity/cns_AUMC_pipeline_v2.def)

Recommended input layout:

- input/<sample>/aligned.bam
- input/<run>/sequencing_summary.txt

Then edit config/samples.yaml to point to your real files.
