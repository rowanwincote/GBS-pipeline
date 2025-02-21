#!/usr/bin/env nextflow

// Define input parameters
params.reads_dir = '~/rw2074/pipeline/reads/'  // Directory with compressed reads
params.kraken_db = '~/rw2074/pipeline/minikraken2_v1_8GB'  // Kraken2 database
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
    .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
    .set { reads_ch }

// Debug: Print detected files
reads_ch.view { sample_id, r1, r2 -> "FOUND: Sample ${sample_id}, Read1: ${r1}, Read2: ${r2}" }

// Kraken2 Process
process kraken2 {
    tag "Kraken2"

    input:
    tuple val(sample_id), path(read1), path(read2)
    val kraken_db  // Use 'val' instead of 'path'

    output:
    path "${sample_id}.kreport", emit: kreport 

    cpus 10

    script:
    """
    echo "Running Kraken2 on ${sample_id}..."
    kraken2 --db ${kraken_db} --paired ${read1} ${read2} --report ${sample_id}.kreport --output /dev/null --memory-mapping --threads 10
    """
}

// Bracken Process
process bracken {
    tag "Bracken"

    input:
    path kreport
    val kraken_db  // Use 'val' instead of 'path'

    output:
    path "${kreport.baseName}.bracken", emit: bracken_output

    script:
    """
    echo "Running Bracken on ${kreport}..."
    bracken -d ${kraken_db} -i ${kreport} -o ${kreport.baseName}.bracken -l ${params.classification_level} -t ${params.threshold}
    """
    bracken_output.view { file -> "Emitting: ${file}" }
}

process generate_summary_report {
    tag "Generate Summary"

    input:
    tuple val(sample_id), path(bracken_output)

    output:
    path "summary_report.txt"

    script:
    """
    echo "Generating summary of GBS abundance..." > summary_report.txt
    for report in *.bracken; do
        if [[ -s "\$report" ]]; then
            echo "-------------------------" >> summary_report.txt
            echo "Sample: \$(basename \$report .bracken)" >> summary_report.txt
            echo "-------------------------" >> summary_report.txt
            cat \$report >> summary_report.txt
            echo "" >> summary_report.txt  # Adds a newline for readability
        else
            echo "WARNING: \$(basename \$report .bracken) has an empty Bracken report!" >> summary_report.txt
        fi
    done
    """
}

// Workflow Execution
workflow {
    kraken_results = kraken2(reads_ch, params.kraken_db)
    bracken_results = bracken(kraken_results.kreport, params.kraken_db)
    generate_summary_report(bracken_results.bracken_output)
}

