"""Maitra et al., 2023"""
import pandas as pd
from os import listdir, rename, getcwd
from os.path import join, basename, dirname, abspath
from pathlib import Path
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("7.20.0")

##### load config and sample sheets #####
configfile: "config.yaml"
samples = pd.read_table(config["samples"]).set_index("Run", drop=False)
resolutions = [0.001]
bprj = "PRJNA883411"
prj  = "maitra2023-dlpfc-mdd"

def plots_doublets_raw(wildcards):
    x = "output/figures/{wildcards.run}_raw/doublets_call".format(wildcards=wildcards)
    return x.replace("\.", "_")


def get_mem_mb(wildcards, attempt):
    return attempt * 500000


##### target rules #####

shell.executable("/bin/bash")

rule all:
    input:
        expand("cellbender/{run}/{run}_output_filtered.h5",
                run=samples["Run"]),
        expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5",
                run=samples["Run"]),
        expand("cellranger/{run}/outs/filtered_feature_bc_matrix.h5",
                run=samples["Run"]),
        expand("scrublet/{run}/{run}_initial_annotation.h5ad",
                run=samples["Run"])

##### load rules #####

CELLRANGER="cd cellranger && source /home/etretiakov/src/cellranger-7.1.0/sourceme.bash && cellranger "

rule cellranger_count:
    input:
        sample=directory("fastq"),
        idx=directory("/home/etretiakov/src/scRNA-seq-references/human_GRCh38_optimized_reference_v2")
    output:
        raw="cellranger/{run}/outs/raw_feature_bc_matrix.h5",
        filtered="cellranger/{run}/outs/filtered_feature_bc_matrix.h5",
        summary="cellranger/{run}/outs/web_summary.html",
        bam="cellranger/{run}/outs/possorted_genome_bam.bam",
    params:
        ids="cellranger/{run}",
        sample="{run}"
    threads: 32
    resources:
        mem_mb=64000
    shell:
        ("{CELLRANGER} count --include-introns true \
            --id={params.sample} \
            --sample={params.sample} \
            --transcriptome={input.idx} \
            --fastqs={input.sample} \
            --jobmode=local \
            --localcores={threads} ")

rule cellbender:
    input:
        "cellranger/{run}/outs/raw_feature_bc_matrix.h5"
    output:
        expand(["cellbender/{{run}}/{{run}}_output.h5", "cellbender/{{run}}/{{run}}_output_filtered.h5"], res=resolutions)
    params:
        ndroplets=lambda wildcards: samples["NTotalDropletsIncluded"][wildcards.run],
        ncells=lambda wildcards: samples["NTotalCells"][wildcards.run],
        h5="cellbender/{run}/{run}_output.h5"
    container:
        "docker://us.gcr.io/broad-dsde-methods/cellbender:latest"
    threads: 4
    resources:
        nvidia_gpu=1,
        mem_mb=10000
    shell:
        ("cellbender remove-background \
            --input {input} \
            --output {params.h5} \
            --cuda \
            --expected-cells {params.ncells} \
            --total-droplets-included {params.ndroplets} \
            --fpr 0.001 \
            --epochs 150")

rule save_h5ad:
    input:
        filt_h5="cellbender/{run}/{run}_output_filtered.h5"
    output:
        dr="cellbender/{run}/{run}_latent_gene_expression.csv",
        h5ad="cellbender/{run}/{run}_filtered.h5ad"
    params:
        sample_run_name="{run}"
    container:
        "docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0"
    threads: 8
    resources:
        mem_mb=20000
    script:
        "code/cb-z.py"

rule doublets_call:
    input:
        h5ad="cellbender/{run}/{run}_filtered.h5ad"
    output:
        scrublet_calls="scrublet/{run}/{run}_scrublet_calls.tsv",
        h5ad="scrublet/{run}/{run}_initial_annotation.h5ad"
    params:
        expected_dblt=lambda wildcards: samples["NExpectedDoubletRate"][wildcards.run],
        plots=plots_doublets_raw
    container:
        "docker://etretiakov/scrna-seq:jammy-2022.12.09-v0.0.1"
    threads: 8
    resources:
        mem_mb=20000
    script:
        "code/scrublet.py"
