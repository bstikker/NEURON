# Portable Snakefile for CNS tumor classification using Nanopore data

configfile: "config/samples.yaml"

import os
from datetime import datetime

SAMPLES = list(config["samples"].keys())
CONTAINER = config.get("container", "singularity/cns_AUMC_pipeline_v2.1.sif")
THREADS = config.get("threads", {})
REFERENCES = config.get("references", {})
TIMECONFIG = config.get("timecourse", {})
TIMECourse_BINS = TIMECONFIG.get("bins", [15, 30, 45, 60, 90, 120, 180, 240, 360, 480])
TIMECourse_CNV_BINS = TIMECONFIG.get("cnv_bins", TIMECourse_BINS)
TIMECourse_CLEANUP = TIMECONFIG.get("cleanup_bins", False)
TIMECourse_CLEANUP_READS = TIMECONFIG.get("cleanup_read_lists", False)

wildcard_constraints:
    sample=r"[^/]+",
    timebin=r"\d+"


def has_timecourse(sample):
    s = config["samples"][sample]
    return all(k in s for k in ["bam", "barcode_label", "sequencing_summary"])


def estimate_mem_mb(wildcards):
    bam_path = config["samples"][wildcards.sample]["bam"]
    if not os.path.exists(bam_path):
        return 32000
    size_gb = os.path.getsize(bam_path) / (1024**3)

    if size_gb < 15:
        mem_gb = 32
    elif size_gb < 30:
        mem_gb = 64
    elif size_gb < 45:
        mem_gb = 96
    else:
        mem_gb = 128

    return mem_gb * 1024


rule all:
    input:
        expand("results/{sample}/merged_probes_methyl_calls_general.csv", sample=SAMPLES),
        expand("results/{sample}/sturgeon_v2_outcome.csv", sample=SAMPLES),
        expand("results/{sample}/sturgeon_v2_outcome.png", sample=SAMPLES),
        expand("results/{sample}/QDNAseq_ACE/{sample}_CNV.png", sample=SAMPLES),
        expand("results/{sample}/QDNAseq_ACE/{sample}_CGHcall_segments.tsv", sample=SAMPLES),
        expand("results/{sample}/QDNAseq_ACE/ACE_summary.tsv", sample=SAMPLES),
        expand("results/{sample}/MGMT_analysis/mgmt_prediction.tsv", sample=SAMPLES),
        expand("results/{sample}/final_report.pdf", sample=SAMPLES)


rule sturgeon:
    input:
        bam=lambda wc: config["samples"][wc.sample]["bam"]
    output:
        modkit_txt=temp("results/{sample}/modkit_extracted.txt"),
        sturgeon_csv="results/{sample}/merged_probes_methyl_calls_general.csv",
        adjusted_bam=temp("results/{sample}/adjusted_merged_sorted.bam"),
        sturgeon_pdf="results/{sample}/merged_probes_methyl_calls_general.pdf",
        bed="results/{sample}/merged_probes_methyl_calls.bed"
    threads: THREADS.get("sturgeon", 10)
    resources:
        mem_mb=estimate_mem_mb,
        runtime=360
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
        bed="results/{sample}/QDNAseq_ACE/{sample}_500kbp.bed",
        seg="results/{sample}/QDNAseq_ACE/{sample}_500kbp.seg",
        plot="results/{sample}/QDNAseq_ACE/{sample}_segmented.png",
        rds="results/{sample}/QDNAseq_ACE/{sample}_copyNumbersSegmented.rds",
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


rule qdnaseq_cghcall_annotate:
    input:
        rds="results/{sample}/QDNAseq_ACE/{sample}_copyNumbersSegmented.rds",
        gene_bed=REFERENCES.get("cnv_gene_bed", "reference/genes/relevant_genes_with_chm13v2_500kb_bin_nrs_fusions_singlebin.bed")
    output:
        cnv="results/{sample}/QDNAseq_ACE/{sample}_CNV.png",
        segments="results/{sample}/QDNAseq_ACE/{sample}_CGHcall_segments.tsv",
        qc="results/{sample}/QDNAseq_ACE/{sample}_QDNAseq_QC.tsv"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=120
    container:
        CONTAINER
    shell:
        r"""
        Rscript scripts/qdnaseq_cghcall_annotate.R \
            {wildcards.sample} \
            {input.rds} \
            {input.gene_bed} \
            results/{wildcards.sample}/QDNAseq_ACE \
            FALSE
        """


rule mgmt_predict_promoter:
    input:
        modkit="results/{sample}/modkit_extracted.txt",
        model=REFERENCES.get("mgmt_model", "scripts/mgmt_137sites_mean_model.Rdata"),
        bed=REFERENCES.get("mgmt_bed", "scripts/mgmt_promoter_coordinates_T2T.bed")
    output:
        result="results/{sample}/MGMT_analysis/mgmt_prediction.tsv"
    resources:
        mem_mb=32000
    container:
        CONTAINER
    shell:
        """
        Rscript scripts/mgmt_analysis.R \
            --modkit {input.modkit} \
            --model {input.model} \
            --bed {input.bed} \
            --out {output.result} \
            --verbose
        """


rule sturgeon_v2:
    input:
        bed="results/{sample}/merged_probes_methyl_calls.bed",
        model="reference/models/cns-v2.zip"
    output:
        csv="results/{sample}/sturgeon_v2_outcome.csv",
        png="results/{sample}/sturgeon_v2_outcome.png"
    params:
        env_bin="/home/P092309/miniconda3/envs/snakemake/bin",
        rscript="scripts/plot_sturgeon_v2.R"
    shell:
        """
        {params.env_bin}/sturgeon-v2 -i {input.bed} -m {input.model} -o {output.csv} -f bed
        Rscript {params.rscript} {output.csv} {wildcards.sample} {output.png}
        """

rule report:
    input:
        sturgeon_pdf="results/{sample}/merged_probes_methyl_calls_general.pdf",
        mgmt="results/{sample}/MGMT_analysis/mgmt_prediction.tsv",
        cnv="results/{sample}/QDNAseq_ACE/{sample}_CNV.png",
        ace="results/{sample}/QDNAseq_ACE/ACE_summary.tsv",
        sturgeon_v2_png="results/{sample}/sturgeon_v2_outcome.png",
        rmd="scripts/report.Rmd"
    output:
        pdf="results/{sample}/final_report.pdf"
    params:
        sample="{sample}",
        version="v1.0.0",
        rundate=lambda wc: datetime.today().strftime("%Y-%m-%d")
    container:
        CONTAINER
    shell:
        r"""
        set -euo pipefail

        OUTDIR="results/{wildcards.sample}"

        rm -f "$OUTDIR"/final_report.aux \
              "$OUTDIR"/final_report.log \
              "$OUTDIR"/final_report.out \
              "$OUTDIR"/final_report.toc

        set +e
        Rscript -e "rmarkdown::render(
            input = '{input.rmd}',
            params = list(
                sample = '{wildcards.sample}',
                version = '{params.version}',
                rundate = '{params.rundate}'
            ),
            output_dir = normalizePath('$OUTDIR', mustWork = TRUE),
            output_file = 'final_report.pdf',
            knit_root_dir = normalizePath('.', mustWork = TRUE)
        )"
        RSTATUS=$?
        set -e

        if [[ $RSTATUS -ne 0 ]]; then
            echo "rmarkdown::render() failed; attempting fallback compile from existing TeX."

            if [[ -f "$OUTDIR/final_report.tex" ]]; then
                (
                    cd "$OUTDIR"
                    xelatex -interaction=nonstopmode -halt-on-error final_report.tex
                    xelatex -interaction=nonstopmode -halt-on-error final_report.tex
                )
            else
                echo "No final_report.tex found, cannot run fallback XeLaTeX compile."
                exit $RSTATUS
            fi
        fi

        test -f "$OUTDIR/final_report.pdf"
        """

rule subset_bam_for_cnv_timepoint:
    input:
        bam=lambda wc: config["samples"][wc.sample]["bam"],
        read_list="read_lists/{sample}/{sample}_read_ids_{timebin}min.txt"
    output:
        bam=temp("results/{sample}/{sample}_{timebin}min/subset_for_cnv.bam"),
        bai=temp("results/{sample}/{sample}_{timebin}min/subset_for_cnv.bam.bai")
    threads: 4
    resources:
        mem_mb=8000,
        runtime=240
    container:
        CONTAINER
    shell:
        r"""
        mkdir -p results/{wildcards.sample}/{wildcards.sample}_{wildcards.timebin}min
        samtools view -@ {threads} -b -N {input.read_list} -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """


rule qdnaseq_ace_timecourse:
    input:
        bam="results/{sample}/{sample}_{timebin}min/subset_for_cnv.bam",
        bai="results/{sample}/{sample}_{timebin}min/subset_for_cnv.bam.bai",
        ace_script="scripts/run_qdnaseq_ace.R",
        cgh_script="scripts/qdnaseq_cghcall_annotate.R",
        gene_bed=REFERENCES.get("cnv_gene_bed", "reference/genes/relevant_genes_with_chm13v2_500kb_bin_nrs_fusions_singlebin.bed")
    output:
        bed="results/{sample}/{sample}_{timebin}min/QDNAseq_ACE/{sample}_{timebin}min_500kbp.bed",
        seg="results/{sample}/{sample}_{timebin}min/QDNAseq_ACE/{sample}_{timebin}min_500kbp.seg",
        cnv_png="results/{sample}/{sample}_{timebin}min/QDNAseq_ACE/{sample}_{timebin}min_CNV.png",
        seg_png="results/{sample}/{sample}_{timebin}min/QDNAseq_ACE/{sample}_{timebin}min_segmented.png",
        rds="results/{sample}/{sample}_{timebin}min/QDNAseq_ACE/{sample}_{timebin}min_copyNumbersSegmented.rds",
        ace_summary="results/{sample}/{sample}_{timebin}min/QDNAseq_ACE/ACE_summary.tsv",
        ace_matrix="results/{sample}/{sample}_{timebin}min/QDNAseq_ACE/{sample}_{timebin}min_ACE_matrixplot.png",
        cgh_segments="results/{sample}/{sample}_{timebin}min/QDNAseq_ACE/{sample}_{timebin}min_CGHcall_segments.tsv",
        qc="results/{sample}/{sample}_{timebin}min/QDNAseq_ACE/{sample}_{timebin}min_QDNAseq_QC.tsv"
    threads: THREADS.get("qdnaseq_ace", 8)
    resources:
        mem_mb=10000,
        runtime=720
    container:
        CONTAINER
    shell:
        r"""
        set -euo pipefail

        OUTDIR="results/{wildcards.sample}/{wildcards.sample}_{wildcards.timebin}min/QDNAseq_ACE"
        mkdir -p "$OUTDIR"

        Rscript {input.ace_script} \
            {wildcards.sample}_{wildcards.timebin}min \
            {input.bam} \
            "$OUTDIR"

        Rscript {input.cgh_script} \
            {wildcards.sample}_{wildcards.timebin}min \
            {output.rds} \
            {input.gene_bed} \
            "$OUTDIR" \
            FALSE
        """


# -----------------------------
# Timecourse workflow
# -----------------------------

SAMPLES = config["samples"].keys()
TIME_BINS = config["timecourse"]["bins"]

rule generate_read_lists:
    input:
        summary=lambda wc: config["samples"][wc.sample]["sequencing_summary"]
    output:
        read_counts="read_lists/{sample}/{sample}_read_counts.tsv",
        read_ids=expand(
            "read_lists/{{sample}}/{{sample}}_read_ids_{timebin}min.txt",
            timebin=config["timecourse"]["bins"]
        )
    params:
        script="scripts/generate_reads_lists.py",
        barcode_label=lambda wc: config["samples"][wc.sample]["barcode_label"],
        bins=lambda wc: " ".join(str(x) for x in config["timecourse"]["bins"])
    container:
        CONTAINER
    shell:
        r"""
        mkdir -p read_lists/{wildcards.sample}
        python {params.script} \
          --summary {input.summary} \
          --sample {wildcards.sample} \
          --barcode-label "{params.barcode_label}" \
          --output read_lists/{wildcards.sample} \
          --bins {params.bins}
        """

rule subset_and_run_sturgeon_timecourse:
    input:
        bam=lambda wc: config["samples"][wc.sample]["bam"],
        read_list="read_lists/{sample}/{sample}_read_ids_{timebin}min.txt",
        probes=REFERENCES.get("sturgeon_probes", "reference/probes/probelocs_chm13.bed"),
        model=REFERENCES.get("sturgeon_model", "reference/models/general.zip")
    output:
        bed="results/{sample}/{sample}_{timebin}min/merged_probes_methyl_calls.bed",
        csv="results/{sample}/{sample}_{timebin}min/merged_probes_methyl_calls_general.csv",
        pdf="results/{sample}/{sample}_{timebin}min/merged_probes_methyl_calls_general.pdf"
    threads: THREADS.get("sturgeon", 10)
    resources:
        mem_mb=estimate_mem_mb,
        runtime=360
    container:
        CONTAINER
    shell:
        """
        mkdir -p results/{wildcards.sample}/{wildcards.sample}_{wildcards.timebin}min
        bash scripts/subset_bam_and_run_sturgeon.sh \
            {input.bam} \
            {input.read_list} \
            {wildcards.sample} \
            {wildcards.timebin} \
            .
        """


rule sturgeon_v2_timecourse:
    input:
        bed="results/{sample}/{sample}_{timebin}min/merged_probes_methyl_calls.bed",
        model="reference/models/cns-v2.zip"
    output:
        csv="results/{sample}/{sample}_{timebin}min/sturgeon_v2_outcome.csv"
    params:
        env_bin="/home/P092309/miniconda3/envs/snakemake/bin"
    shell:
        """
        {params.env_bin}/sturgeon-v2 -i {input.bed} -m {input.model} -o {output.csv} -f bed
        """

rule timecourse_rmd:
    input:
        sturgeon_v1_csvs=expand(
            "results/{{sample}}/{{sample}}_{timebin}min/merged_probes_methyl_calls_general.csv",
            timebin=TIME_BINS
        ),
        sturgeon_v2_csvs=expand(
            "results/{{sample}}/{{sample}}_{timebin}min/sturgeon_v2_outcome.csv",
            timebin=TIME_BINS
        ),
        cnv_pngs=expand(
            "results/{{sample}}/{{sample}}_{timebin}min/QDNAseq_ACE/{{sample}}_{timebin}min_CNV.png",
            timebin=TIME_BINS
        ),
        read_counts="read_lists/{sample}/{sample}_read_counts.tsv"
    output:
        pdf="results/{sample}/{sample}_timecourse.pdf"
    params:
        sample_id="{sample}"
    container:
        CONTAINER
    shell:
        r"""
        Rscript -e "rmarkdown::render(
            input = 'scripts/timecourse_report.Rmd',
            params = list(
                sample_id = '{wildcards.sample}'
            ),
            output_dir = normalizePath('results/{wildcards.sample}', mustWork = TRUE),
            output_file = '{wildcards.sample}_timecourse.pdf',
            knit_root_dir = normalizePath('.', mustWork = TRUE)
        )"
        """

rule all_timecourse:
    input:
        pdfs=expand("results/{sample}/{sample}_timecourse.pdf", sample=SAMPLES)
