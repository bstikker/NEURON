# Portable Snakefile for CNS tumor classification using Nanopore data

configfile: "config/samples.yaml"

SAMPLES = list(config["samples"].keys())
CONTAINER = config.get("container", "singularity/cns_AUMC_pipeline_v2.sif")
THREADS = config.get("threads", {})
REFERENCES = config.get("references", {})

rule all:
    input:
        expand("results/{sample}/merged_probes_methyl_calls_general.csv", sample=SAMPLES),
        expand("results/{sample}/QDNAseq_ACE/{sample}_segmented.png", sample=SAMPLES),
        expand("results/{sample}/QDNAseq_ACE/ACE_summary.tsv", sample=SAMPLES),
        expand("results/{sample}/MGMT_analysis/mgmt_prediction.tsv", sample=SAMPLES),
        expand("results/{sample}/final_report.pdf", sample=SAMPLES)

rule sturgeon:
    input:
        bam=lambda wildcards: config["samples"][wildcards.sample]["bam"]
    output:
        modkit_txt=temp("results/{sample}/modkit_extracted.txt"),
        sturgeon_csv="results/{sample}/merged_probes_methyl_calls_general.csv",
        adjusted_bam=temp("results/{sample}/adjusted_merged_sorted.bam"),
        sturgeon_pdf="results/{sample}/merged_probes_methyl_calls_general.pdf"
    resources:
        mem_mb=20000,
        threads=THREADS.get("sturgeon", 10)
    params:
        sample="{sample}"
    container:
        CONTAINER
    shell:
        """
        Rscript scripts/run_sturgeon.R {input.bam} {params.sample}
        """

rule qdnaseq_ace:
    input:
        bam="results/{sample}/adjusted_merged_sorted.bam"
    output:
        bed=temp("results/{sample}/QDNAseq_ACE/{sample}_500kbp.bed"),
        seg=temp("results/{sample}/QDNAseq_ACE/{sample}_500kbp.seg"),
        plot="results/{sample}/QDNAseq_ACE/{sample}_segmented.png",
        rds=temp("results/{sample}/QDNAseq_ACE/{sample}_copyNumbersSegmented.rds"),
        ace_summary="results/{sample}/QDNAseq_ACE/ACE_summary.tsv",
        ace_plot="results/{sample}/QDNAseq_ACE/{sample}_ACE_matrixplot.png"
    params:
        sample="{sample}"
    resources:
        mem_mb=10000,
        threads=THREADS.get("qdnaseq_ace", 8)
    container:
        CONTAINER
    shell:
        """
        Rscript scripts/run_qdnaseq_ace.R {params.sample} {input.bam}
        """

rule mgmt_predict_promoter:
    input:
        modkit="results/{sample}/modkit_extracted.txt",
        model=REFERENCES.get("mgmt_model", "scripts/mgmt_137sites_mean_model.Rdata"),
        bed=REFERENCES.get("mgmt_bed", "scripts/mgmt_promoter_coordinates_T2T.bed")
    output:
        result="results/{sample}/MGMT_analysis/mgmt_prediction.tsv"
    container:
        CONTAINER
    shell:
        """
        Rscript scripts/mgmt_analysis.R             --modkit {input.modkit}             --model {input.model}             --bed {input.bed}             --out {output.result}             --verbose
        """

rule report:
    input:
        sturgeon_pdf="results/{sample}/merged_probes_methyl_calls_general.pdf",
        mgmt="results/{sample}/MGMT_analysis/mgmt_prediction.tsv",
        cnv_png="results/{sample}/QDNAseq_ACE/{sample}_segmented.png",
        ace="results/{sample}/QDNAseq_ACE/ACE_summary.tsv",
        rmd="scripts/report.Rmd"        
    output:
        pdf="results/{sample}/final_report.pdf"
    params:
        outdir=lambda wildcards: f"results/{wildcards.sample}",
        sample="{sample}"
    container:
        CONTAINER
    shell:
        r"""
        Rscript -e "rmarkdown::render(\
          input = 'scripts/report.Rmd', \
          params = list(sample = '{wildcards.sample}'), \
          output_dir = normalizePath('results/{wildcards.sample}', mustWork = TRUE), \
          output_file = 'final_report.pdf', \
          knit_root_dir = normalizePath('.', mustWork = TRUE) \
        )"
        """
