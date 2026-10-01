# Reference genome

The pipeline does not read a reference genome: it starts from BAM files that are already aligned, and no step needs the FASTA.

Input BAM files must be aligned to **T2T-CHM13 v2.0**, because the Sturgeon probe locations, QDNAseq bins, MGMT promoter coordinates and gene annotations all use T2T-CHM13 v2.0 coordinates.

FASTA download, if you need it for alignment: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
