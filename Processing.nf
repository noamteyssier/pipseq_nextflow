#!/usr/bin/env nextflow

workflow {
  read_pairs = Channel
    .fromFilePairs ( params.reads, checkIfExists: true )

  PIPSpeak(
    read_pairs
  )

  KallistoBUS(
    PIPSpeak.out.filtered_ch
  )

  BUStoolsSort( 
    KallistoBUS.out.bus_ch
  )

  BUStoolsInspect(
    BUStoolsSort.out.sorted_bus_ch,
    KallistoBUS.out.ec_ch
  )

  BUStoolsCount( 
    BUStoolsSort.out.sorted_bus_ch,
    KallistoBUS.out.transcripts_ch,
    KallistoBUS.out.ec_ch
  )

  BuildH5AD(
    BUStoolsCount.out.mtx_ch,
    BUStoolsCount.out.barcodes_ch,
    BUStoolsCount.out.genes_ch
  )

  PlotQC(
    BuildH5AD.out.h5ad_ch
  )

  FilterNormLog(
    BuildH5AD.out.h5ad_ch
  )

  Merge(
    FilterNormLog.out.h5ad_filt_ch.collect(),
  )

  SingleCell(
    Merge.out.h5ad_ch
  )
}

process PIPSpeak {

    publishDir "${params.outdir}/pipspeak/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("filtered_${sample_id}_R[12].fq.gz"), emit: filtered_ch
    path("filtered_${sample_id}_whitelist.txt")
    path("filtered_${sample_id}_log.yaml")

    script:
    """
    #!/bin/bash

    ## touch filtered_${sample_id}_R1.fq.gz
    ## touch filtered_${sample_id}_R2.fq.gz
    ## touch filtered_${sample_id}_whitelist.txt
    ## touch filtered_${sample_id}_log.yaml

    pipspeak \
        -c ${params.pipspeak_yaml} \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -p filtered_${sample_id} \
        -t ${params.n_threads}
    """
}

process KallistoBUS {
    
    publishDir "${params.outdir}/kallisto_bus/", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}/output.bus"), emit: bus_ch
    path("${sample_id}/matrix.ec"), emit: ec_ch
    path("${sample_id}/transcripts.txt"), emit: transcripts_ch

    script:
    """
    #!/bin/bash

    ## mkdir -p ${sample_id}
    ## touch ${sample_id}/output.bus
    ## touch ${sample_id}/matrix.ec
    ## touch ${sample_id}/transcripts.txt

    kallisto bus -n -t ${params.n_threads} \
        -i ${params.kallisto_index} \
        -o ${sample_id} \
        -x ${params.pipseq_tech} \
        ${reads[0]} ${reads[1]}
    """
}

process BUStoolsSort {
    
    publishDir "${params.outdir}/kallisto_bus/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(unsorted_bus)

    output:
    tuple val(sample_id), path("output.sorted.bus"), emit: sorted_bus_ch

    script:
    """
    #!/bin/bash

    ## mkdir -p ${sample_id}
    ## touch ${sample_id}/output.sorted.bus

    bustools sort -t ${params.n_threads} \
        -o output.sorted.bus \
        ${unsorted_bus}
    """
}

process BUStoolsInspect {
    
    publishDir "${params.outdir}/kallisto_bus/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(sorted_bus)
    path(ec)

    output:
    path("inspect.json")

    script:
    """
    #!/bin/bash

    ## touch inspect.json

    bustools inspect \
        -o inspect.json \
        -e ${ec} \
        ${sorted_bus}
    """

}

process BUStoolsCount {

    publishDir "${params.outdir}/counts/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(sorted_bus)
    path(transcripts)
    path(ec)

    output:
    tuple val(sample_id), path("data.mtx"), emit: mtx_ch
    path "data.barcodes.txt", emit: barcodes_ch
    path "data.genes.txt", emit: genes_ch

    script:
    """
    #!/bin/bash

    ## mkdir -p ${sample_id}
    ## touch ${sample_id}/data.mtx
    ## touch ${sample_id}/data.barcodes.txt
    ## touch ${sample_id}/data.genes.txt

    bustools count \
        -o data \
        -t ${transcripts} \
        -g ${params.t2g} \
        -e ${ec} \
        --genecounts \
        ${sorted_bus}
    """
}

process BuildH5AD {
    
    publishDir "${params.outdir}/counts/${sample_id}", mode: 'symlink'
    conda "${params.scanpy_env}"

    input:
    tuple val(sample_id), path(mtx)
    path(barcodes)
    path(genes)

    output:
    tuple val(sample_id), path("${sample_id}.h5ad"), emit: h5ad_ch

    script:
    """
    #!/usr/bin/env python3

    import scanpy as sc
    import anndata as ad

    # Load data
    adata = sc.read_mtx("${mtx}")
    barcodes = open("${barcodes}").read().splitlines()
    genes = open("${genes}").read().splitlines()

    # Set index
    adata.obs.index = barcodes
    adata.var.index = genes

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Add SampleID to obs
    adata.obs["sample_id"] = "${sample_id}"

    # Write to file
    adata.write_h5ad("${sample_id}.h5ad")
    """
}

process PlotQC {
    
    publishDir "${params.outdir}/qc/${sample_id}", mode: 'symlink'
    conda "${params.scanpy_env}"

    input:
    tuple val(sample_id), path(h5ad)

    output:
    path("rankplot.svg")

    script:
    """
    #!/usr/bin/env python3

    import numpy as np
    import scanpy as sc
    import matplotlib.pyplot as plt

    # Load data
    adata = sc.read_h5ad("${h5ad}")

    # Generate Axis
    knee = np.sort(adata.X.toarray().sum(axis=1).flatten())[::-1]
    fig, ax = plt.subplots(figsize=(10, 7))

    ax.loglog(knee, range(len(knee)),linewidth=5, color="g")
    ax.axvline(${params.cell_umi_threshold}, color="r", linestyle="--", linewidth=3)

    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Set of Barcodes")

    plt.grid(True, which="both")
    plt.tight_layout()
    plt.savefig("rankplot.svg")
    """
}

process FilterNormLog {
    
    publishDir "${params.outdir}/counts/${sample_id}", mode: 'symlink'
    conda "${params.scanpy_env}"

    input:
    tuple val(sample_id), path(h5ad)

    output:
    path("${sample_id}.filt.h5ad"), emit: h5ad_filt_ch
    val(sample_id), emit: sample_id_ch

    script:
    """
    #!/usr/bin/env python3

    import scanpy as sc

    # Load data
    adata = sc.read_h5ad("${h5ad}")

    # Filter cells
    sc.pp.filter_cells(adata, min_counts=${params.cell_umi_threshold})

    # Filter genes
    sc.pp.filter_genes(adata, min_counts=${params.gene_cell_threshold})

    # Set Raw Data
    adata.raw = adata

    # Normalize and Log
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Write to file
    adata.write_h5ad("${sample_id}.filt.h5ad")
    """
}

process Merge {

    publishDir "${params.outdir}/counts", mode: 'symlink'
    conda "${params.scanpy_env}"

    input:
    path(adata_filt_ch)

    output:
    path("adata.h5ad"), emit: h5ad_ch

    script:
    """
    #!/usr/bin/env python3

    import scanpy as sc

    def build_g2s(t2g_path: str) -> dict:
        used = dict()
        g2s = dict()
        for line in open(t2g_path):
            record = line.strip().split('\t')
            tx = record[1]
            sym = record[2]
            if tx not in g2s:
                if sym == "":
                    g2s[tx] = tx
                else:
                    if sym not in used:
                        used[sym] = 0
                    else:
                        used[sym] += 1
                        sym += f"_{used[sym]}"
                    g2s[tx] = sym
        return g2s

    def concatenate_h5ads(path_list: list):
        adata = sc.concat(
            [sc.read_h5ad(i) for i in path_list],
            join='outer',
        )
        adata.obs_names_make_unique()
        return adata


    path_list = [i for i in "${adata_filt_ch}".split(" ")]
    adata = concatenate_h5ads(path_list)
    g2s = build_g2s("${params.t2g}")
    
    adata.var["ensembl_id"] = adata.var.index
    adata.var["gene_symbol"] = adata.var["ensembl_id"].map(g2s)
    adata.var.index = adata.var["gene_symbol"]
    adata.var_names_make_unique()

    adata.write_h5ad("adata.h5ad")
    """
}

process SingleCell {

    publishDir "${params.outdir}/analysis", mode: 'symlink'
    conda "${params.scanpy_env}"

    input:
    path(merged_h5ad)

    output:
    path(merged_h5ad), emit: sc_h5ad_ch

    script:
    """
    #!/usr/bin/env python3

    import scanpy as sc

    # Load data
    adata = sc.read_h5ad("${merged_h5ad}")

    # Calculate Highly Variable Genes
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    # Perform PCA
    sc.tl.pca(adata)

    # Perform Nearest Neighbor Graph
    sc.pp.neighbors(adata)

    # Perform UMAP
    sc.tl.umap(adata)

    # Perform Leiden Clustering
    sc.tl.leiden(adata)

    # Write to file
    adata.write_h5ad("${merged_h5ad}")
    """

}
