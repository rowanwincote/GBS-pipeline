#!/usr/bin/env nextflow

// Define input parameters
params.read1 = '/mnt/c/my_nextflow_work/pipeline/20280_5#33_1.fastq'
params.read2 = '/mnt/c/my_nextflow_work/pipeline/20280_5#33_2.fastq'
params.kraken_db = '/mnt/c/my_nextflow_work/pipeline/minikraken2_v1_8GB'
params.bracken_db = '/mnt/c/my_nextflow_work/pipeline/Bracken'
params.classification_level = 'S'
params.threshold = 10
params.read_len = 150

process kraken2 {
    tag "Kraken2"
    input:
    tuple path(read1), path(read2)

    output:
    path "${read1.simpleName}.kreport", emit: kreport

    script:
    """
    echo "Running Kraken2 on ${read1} and ${read2}..."
    ls -lh ${read1} ${read2}  # Debugging: Check file existence
    kraken2 --db ${params.kraken_db} --paired ${read1} ${read2} --report ${read1.simpleName}.kreport --output /dev/null
    ls -lh ${read1.simpleName}.kreport  # Debugging: Check if output was created
    """
}

process bracken {
    tag "Bracken"
    input:
    path kreport

    output:
    path "${kreport.baseName}.bracken", emit: bracken_output

    script:
    """
    echo "Running Bracken on ${kreport}..."
    python /mnt/c/my_nextflow_work/pipeline/Bracken/src/est_abundance.py -i ${kreport} \
        -k ${params.kraken_db}/database${params.read_len}mers.kmer_distrib \
        -o ${kreport.baseName}.bracken \
        -l ${params.classification_level} -t ${params.threshold}
    """
}

process generate_summary_report {
    tag "Generate Summary"
    input:
    path bracken_reports

    output:
    path "summary_report.txt"

    script:
    """
    echo "Generating summary of GBS abundance..." > summary_report.txt
    for report in \$(ls ${bracken_reports}); do
        echo "Processing \${report}..." >> summary_report.txt
        cat \${report} >> summary_report.txt
    done
    """
}

workflow {
    if (!file(params.read1).exists() || !file(params.read2).exists()) {
        throw new Exception("ERROR: One or both input files do not exist.")
    }
    if (!file(params.kraken_db).exists()) {
        throw new Exception("ERROR: Kraken2 database not found.")
    }
    if (!file(params.bracken_db).exists()) {
        throw new Exception("ERROR: Bracken database not found.")
    }

    reads = tuple(file(params.read1), file(params.read2))

    kraken_results = kraken2(reads)
    bracken_results = bracken(kraken_results)  // Pass full output
    generate_summary_report(bracken_results)  // Use full channel
}

