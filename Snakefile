# Portable Snakefile for CNS tumor classification using Nanopore data

configfile: "config/samples.yaml"

SAMPLES = list(config["samples"].keys())
CONTAINER = config.get("container", "singularity/cns_AUMC_pipeline_v2.sif")
THREADS = config.get("threads", {})
REFERENCES = config.get("references", {})

import os

TIMECONFIG = config.get("timecourse", {})
TIMECourse_BINS = TIMECONFIG.get("bins", [15, 30, 45, 60, 90, 120, 180, 240, 360, 480])
TIMECourse_CNV_BINS = TIMECONFIG.get("cnv_bins", TIMECourse_BINS)
TIMECourse_CLEANUP = TIMECONFIG.get("cleanup_bins", False)
TIMECourse_CLEANUP_READS = TIMECONFIG.get("cleanup_read_lists", False)

def has_timecourse(sample):
    s = config["samples"][sample]
    return all(k in s for k in ["bam", "barcode_label", "sequencing_summary"])

def estimate_mem_mb(wildcards):
    bam_path = config["samples"][wildcards.sample]["bam"]
    if not os.path.exists(bam_path):
        return 32000
    size_gb = os.path.getsize(bam_path) / (1024**3)

    # more conservative for large BAMs
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
        expand("results/{sample}/QDNAseq_ACE/{sample}_CNV.png", sample=SAMPLES),
        expand("results/{sample}/QDNAseq_ACE/{sample}_CGHcall_segments.tsv", sample=SAMPLES),
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
    threads: THREADS.get("sturgeon", 10)
    resources:
        mem_mb=estimate_mem_mb,  # dynamically estimate per sample
        runtime=360  # keep 6h max
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
        gene_bed="reference/genes/relevant_genes_with_chm13v2_500kb_bin_nrs_fusions_singlebin.bed"
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
        Rscript scripts/mgmt_analysis.R             --modkit {input.modkit}             --model {input.model}             --bed {input.bed}             --out {output.result}             --verbose
        """

rule report:
    input:
        sturgeon_pdf="results/{sample}/merged_probes_methyl_calls_general.pdf",
        mgmt="results/{sample}/MGMT_analysis/mgmt_prediction.tsv",
        cnv="results/{sample}/QDNAseq_ACE/{sample}_CNV.png",
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

rule all_timecourse:
    input:
        expand("results/{sample}/{sample}_timecourse_report.pdf",
               sample=[s for s in SAMPLES if has_timecourse(s)]),
        expand("results/{sample}/{sample}_timecourse_summary.tsv",
               sample=[s for s in SAMPLES if has_timecourse(s)]),
        expand("read_lists/{sample}/{sample}_read_counts.tsv",
               sample=[s for s in SAMPLES if has_timecourse(s)]),

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
        gene_bed="reference/genes/relevant_genes_with_chm13v2_500kb_bin_nrs_fusions_singlebin.bed"
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

rule sturgeon_timecourse:
    input:
        bam=lambda wc: config["samples"][wc.sample]["bam"],
        summary=lambda wc: config["samples"][wc.sample]["sequencing_summary"],
        script_py="scripts/generate_reads_lists.py",
        script_sh="scripts/subset_bam_and_run_sturgeon.sh",
        script_r="scripts/aggregate_sturgeon_timecourse_line_top10.R",
        probes="reference/probes/probelocs_chm13.bed",
        model="reference/models/general.zip",
        cnv_pngs=expand(
            "results/{sample}/{sample}_{timebin}min/QDNAseq_ACE/{sample}_{timebin}min_CNV.png",
            sample=lambda wc: wc.sample,
            timebin=TIMECourse_CNV_BINS
        )
    output:
        pdf="results/{sample}/{sample}_timecourse_report.pdf",
        tsv="results/{sample}/{sample}_timecourse_summary.tsv",
        read_counts="read_lists/{sample}/{sample}_read_counts.tsv"
    params:
        barcode_label=lambda wc: config["samples"][wc.sample]["barcode_label"],
        bins=lambda wc: " ".join(map(str, TIMECourse_BINS)),
        project_root=lambda wc: workflow.basedir,
        cleanup_bins=lambda wc: "true" if TIMECourse_CLEANUP else "false"
    threads: THREADS.get("sturgeon", 10)
    resources:
        mem_mb=estimate_mem_mb,
        runtime=720
    container:
        CONTAINER
    shell:
        r"""
        mkdir -p read_lists/{wildcards.sample} results/{wildcards.sample} tmp

        python3 {input.script_py} \
            --summary {input.summary} \
            --sample {wildcards.sample} \
            --barcode-label "{params.barcode_label}" \
            --output read_lists/{wildcards.sample} \
            --bins {params.bins}

        for BIN in {params.bins}; do
            READ_LIST="read_lists/{wildcards.sample}/{wildcards.sample}_read_ids_${{BIN}}min.txt"
            if [[ -s "${{READ_LIST}}" ]]; then
                bash {input.script_sh} \
                    "{input.bam}" \
                    "${{READ_LIST}}" \
                    "{wildcards.sample}" \
                    "${{BIN}}" \
                    "{params.project_root}"
            fi
        done

        Rscript {input.script_r} \
            results/{wildcards.sample} \
            {output.pdf} \
            {output.tsv} \
            {output.read_counts}

        if [[ "{params.cleanup_bins}" == "true" ]]; then
            find results/{wildcards.sample} -maxdepth 1 -type d -name "{wildcards.sample}_*min" -exec rm -rf {{}} +
        fi
        """
