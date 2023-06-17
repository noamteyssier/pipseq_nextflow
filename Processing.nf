#!/usr/bin/env nextflow

params.pipspeak_yaml = "data/config_v3.yaml"
params.sequence_dir = "data/sequences"
params.outdir = "results"
params.reads = "$params.sequence_dir/*_R{1,2}*.fastq.gz"
params.n_threads = 2

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

    touch filtered_${sample_id}_R1.fq.gz
    touch filtered_${sample_id}_R2.fq.gz
    touch filtered_${sample_id}_whitelist.txt
    touch filtered_${sample_id}_log.yaml

    # pipspeak \
    #     -c ${params.pipspeak_yaml} \
    #     -i ${reads[0]} \
    #     -I ${reads[1]} \
    #     -p filtered_${sample_id} \
    #     -t ${params.n_threads}
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

    mkdir -p ${sample_id}
    touch ${sample_id}/output.bus
    touch ${sample_id}/matrix.ec
    touch ${sample_id}/transcripts.txt
    """
}

process BUStoolsSort {
    
    publishDir "${params.outdir}/kallisto_bus/", mode: 'symlink'

    input:
    tuple val(sample_id), path(unsorted_bus)

    output:
    tuple val(sample_id), path("${sample_id}/output.sorted.bus"), emit: sorted_bus_ch

    script:
    """
    #!/bin/bash

    mkdir -p ${sample_id}
    touch ${sample_id}/output.sorted.bus
    """
}

process BUStoolsCount {

    publishDir "${params.outdir}/counts/", mode: 'symlink'

    input:
    tuple val(sample_id), path(sorted_bus)
    path(transcripts)
    path(ec)

    output:
    tuple val(sample_id), path("${sample_id}/data.mtx"), emit: mtx_ch
    path "${sample_id}/data.barcodes.txt"
    path "${sample_id}/data.ec.txt"

    script:
    """
    #!/bin/bash

    mkdir -p ${sample_id}
    touch ${sample_id}/data.mtx
    touch ${sample_id}/data.barcodes.txt
    touch ${sample_id}/data.ec.txt
    """

}

workflow {
  read_pairs = Channel
    .fromFilePairs ( params.reads, checkIfExists: true )

  read_pairs | PIPSpeak
  PIPSpeak.out.filtered_ch | KallistoBUS
  KallistoBUS.out.bus_ch | BUStoolsSort
  BUStoolsCount( 
    BUStoolsSort.out.sorted_bus_ch,
    KallistoBUS.out.transcripts_ch,
    KallistoBUS.out.ec_ch
  )
}
