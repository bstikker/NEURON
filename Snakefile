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
        expand("results/{sample}/pseudotime/confidence_vs_pseudotime.png", sample=SAMPLES),
        expand("results/{sample}/pseudotime/confidence_vs_pseudotime.csv", sample=SAMPLES),
        expand("results/{sample}/MGMT_analysis/mgmt_prediction.tsv", sample=SAMPLES),
        expand("results/{sample}/pseudotime/cleanup_done.txt", sample=SAMPLES),
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

rule pseudotime:
    input:
        bam=lambda wildcards: config["samples"][wildcards.sample]["bam"],
        seq_summary=lambda wildcards: config["samples"][wildcards.sample]["summary"],
        modkit="results/{sample}/modkit_extracted.txt"
    output:
        confidence_csv="results/{sample}/pseudotime/confidence_vs_pseudotime.csv",
        confidence_plot="results/{sample}/pseudotime/confidence_vs_pseudotime.png"
    resources:
        mem_mb=32000,
        threads=THREADS.get("pseudotime", 10),
        cpus=THREADS.get("pseudotime", 10)
    params:
        sample="{sample}"
    container:
        CONTAINER
    shell:
        """
        Rscript scripts/pseudotime_analysis.R {input.bam} {params.sample} {input.seq_summary}
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

rule cleanup_bins:
    input:
        "results/{sample}/pseudotime/confidence_vs_pseudotime.png"
    output:
        touch("results/{sample}/pseudotime/cleanup_done.txt")
    shell:
        """
        rm -rf results/{wildcards.sample}/pseudotime/bin_*
        touch {output}
        """

rule estimate_confident_time:
    input:
        csv="results/{sample}/pseudotime/confidence_vs_pseudotime.csv"
    output:
        "results/{sample}/pseudotime/confident_time.csv"
    params:
        label="{sample}"
    resources:
        mem_mb=2000,
        cpus=1
    container:
        CONTAINER
    shell:
        r"""
        python scripts/estimate_confident_time.py           --in_csv {input.csv}           --out_csv {output}           --label {params.label}
        """

rule estimate_stable_time:
    input:
        csv="results/{sample}/pseudotime/confidence_vs_pseudotime.csv"
    output:
        "results/{sample}/pseudotime/stable_time.csv"
    params:
        label="{sample}"
    resources:
        mem_mb=2000,
        cpus=1
    container:
        CONTAINER
    shell:
        r"""
        python scripts/estimate_stable_time.py           --in_csv {input.csv}           --out_csv {output}           --label {params.label}
        """

rule report:
    input:
        sturgeon_pdf="results/{sample}/merged_probes_methyl_calls_general.pdf",
        mgmt="results/{sample}/MGMT_analysis/mgmt_prediction.tsv",
        cnv_png="results/{sample}/QDNAseq_ACE/{sample}_segmented.png",
        ace="results/{sample}/QDNAseq_ACE/ACE_summary.tsv",
        pseudo="results/{sample}/pseudotime/confidence_vs_pseudotime.png",
        confident="results/{sample}/pseudotime/confident_time.csv",
        stable="results/{sample}/pseudotime/stable_time.csv",
        rmd="scripts/report.Rmd",
        logo=REFERENCES.get("report_logo", "scripts/logo.png"),
        latex_header=REFERENCES.get("report_header", "scripts/logo.tex")
    output:
        pdf="results/{sample}/final_report.pdf"
    params:
        outdir=lambda wildcards: f"results/{wildcards.sample}",
        sample="{sample}"
    container:
        CONTAINER
    shell:
        r"""
        cp {input.logo} {params.outdir}/logo.png
        Rscript -e "rmarkdown::render(          input = '{input.rmd}',           params = list(sample = '{params.sample}'),           output_dir = normalizePath('{params.outdir}', mustWork = TRUE),           output_file = 'final_report.pdf',           knit_root_dir = normalizePath('.', mustWork = TRUE),           output_options = list(includes = list(in_header = normalizePath('{input.latex_header}', mustWork = TRUE)))        )"
        """
