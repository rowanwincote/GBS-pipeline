#!/usr/bin/env nextflow

// Define input parameters (merged from main.nf and shovel.nf)
params.read1 = '/mnt/c/my_nextflow_work/pipeline/reads/SRR8541972_1.fastq.gz' // Replace with paths to your FASTQ files
params.read2 = '/mnt/c/my_nextflow_work/pipeline/reads/SRR8541972_2.fastq.gz'
params.kraken_db = '/mnt/c/my_nextflow_work/pipeline/minikraken2_v1_8GB'
params.bracken_db = '/mnt/c/my_nextflow_work/pipeline/Bracken' // Base path for Bracken
params.classification_level = 'S'
params.threshold = 10
params.read_len = 150
params.threads = 8 

process kraken2 {
    tag "Kraken2 for ${read1.simpleName}"
    input:
    tuple path(read1), path(read2)
    output:
    path "${read1.simpleName}.kreport", emit: kreport
    script:
    """
    echo "Running Kraken2 on ${read1} and ${read2}..."
    ls -lh ${read1} ${read2}  # Debugging: Check file existence
    kraken2 --db ${params.kraken_db} --paired ${read1} ${read2} --report ${read1.simpleName}.kreport --output /dev/null --threads ${params.threads}
    ls -lh ${read1.simpleName}.kreport  # Debugging: Check if output was created
    """
}

process bracken {
    tag "Bracken for ${kreport.baseName}"
    input:
    path kreport
    output:
    path "${kreport.baseName}.bracken", emit: bracken_output
    script:

    """
    echo "Running Bracken on ${kreport}..."
    python ${params.bracken_db}/src/est_abundance.py -i ${kreport} \
        -k ${params.kraken_db}/database${params.read_len}mers.kmer_distrib \
        -o ${kreport.baseName}.bracken \
        -l ${params.classification_level} -t ${params.threshold}
    """
}

process generate_summary_report {
    tag "Generate Summary Report"
    input:
    path bracken_report_file 
    output:
    path "summary_report.txt"
    script:
    """
    echo "Generating summary of GBS abundance..." > summary_report.txt
    echo "Processing ${bracken_report_file}..." >> summary_report.txt
    cat ${bracken_report_file} >> summary_report.txt
    """
}

// Processes from shovel.nf
process shovill {
    tag "Shovill for ${read1.simpleName}"
    input:
    tuple path(read1), path(read2)

    output:
    path "shovill_output/contigs.fa", emit: contigs
    path "shovill_output", emit: shovill_dir 

    script:
    """
    shovill --R1 '${read1}' --R2 '${read2}' --outdir shovill_output --cpus ${params.threads}
    """
}

process quast {
    tag "QUAST for assembled contigs"
    input:
    path contigs_file 

    output:
    path "quast_output"

    script:
    """
    quast.py '${contigs_file}' -o quast_output --threads ${params.threads}
    """
}

// Combined Workflow
workflow {
    // Input validation (from main.nf)
    if (!file(params.read1).exists() || !file(params.read2).exists()) {
        error("ERROR: One or both input read files do not exist: ${params.read1}, ${params.read2}")
    }
    if (!file(params.kraken_db).exists()) {
        error("ERROR: Kraken2 database not found: ${params.kraken_db}")
    }
    if (!file(params.bracken_db).exists()) { 
        error("ERROR: Bracken base directory not found: ${params.bracken_db}")
    }
   


    read_pair_input = tuple(file(params.read1), file(params.read2))

    // Path 1: Kraken2 -> Bracken -> Summary (from main.nf)
    kreport_ch = kraken2(read_pair_input)
    bracken_output_ch = bracken(kreport_ch.kreport)
    summary_report_ch = generate_summary_report(bracken_output_ch.bracken_output)

    // Path 2: Shovill -> Quast (from shovel.nf)
    shovill_results_ch = shovill(read_pair_input)
    quast_results_ch = quast(shovill_results_ch.contigs)


}
