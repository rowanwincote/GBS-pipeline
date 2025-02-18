#!/usr/bin/env nextflow

// Define input parameters
params.reads_dir = '/mnt/c/my_nextflow_work/pipeline/reads/'  // Directory with compressed reads
params.kraken_db = '/mnt/c/my_nextflow_work/pipeline/minikraken2_v1_8GB'
params.classification_level = 'S'
params.threshold = 10
params.read_len = 150

// Debug: Check if directory exists
if (!file(params.reads_dir).exists()) {
    error "ERROR: Reads directory '${params.reads_dir}' does not exist!"
}

// Create a channel with all FASTQ.GZ file pairs
Channel
    .fromFilePairs("${params.reads_dir}/*_{1,2}.fastq.gz", checkIfExists: true)
    .ifEmpty { error "ERROR: No FASTQ.GZ files found in '${params.reads_dir}'!" }
    .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }  // FIX: Unpack the list
    .set { reads_ch }

// Debug: Print detected files
reads_ch.view { sample_id, r1, r2 -> "FOUND: Sample ${sample_id}, Read1: ${r1}, Read2: ${r2}" }

// Kraken2 Process
process kraken2 {
    tag "Kraken2"
    maxForks 1
    input:
    tuple val(sample_id), path(read1), path(read2)
    output:
    path "${sample_id}.kreport", emit: kreport 
    cpus 10         
    script:
    """
    echo "Running Kraken2 on ${sample_id}..."
    kraken2 --db ${params.kraken_db} --paired ${read1} ${read2} --report ${sample_id}.kreport --output /dev/null --memory-mapping --threads 10
    """
}

// Bracken Process
process bracken {
    tag "Bracken"

    input:
    path kreport

    output:
    path "${kreport.baseName}.bracken", emit: bracken_output

    script:
    """
    echo "Running Bracken on ${kreport}..."
    bracken -d ${params.kraken_db} -i ${kreport} -o ${kreport.baseName}.bracken -l ${params.classification_level} -t ${params.threshold}
    """
}

// Generate Summary Report
process generate_summary_report {
    tag "Generate Summary"

    input:
    path bracken_output

    output:
    path "summary_report.txt"

    script:
    """
    echo "Generating summary of GBS abundance..." > summary_report.txt
    for report in ${bracken_output}; do
        echo "Processing ${report}..." >> summary_report.txt
        cat ${report} >> summary_report.txt
    done
    """
}

// Workflow Execution
workflow {
    kraken_results = kraken2(reads_ch)
    bracken_results = bracken(kraken_results.kreport)
    generate_summary_report(bracken_results.bracken_output)
}
