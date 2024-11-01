process PREPROCESS {
    tag "$meta.id"
    label 'process_fixed_medium_cpu'
    label 'publish_outdir'

    input:
    tuple val(meta), path(fastq)
    path(fasta)

    output:
    tuple val(meta), path(fastq_R1), path(fastq_R2), emit: fastq
    tuple path(fasta_bgzip), path(fasta_index), emit: fasta

    script:
    def sample_name = meta.name

    fasta_raw   = "fasta.raw"
    fasta_bgzip = "fasta.bgz"
    fasta_index = "${fasta_bgzip}.fai"

    fastq_R1 = "${sample_name}_R1.fastq.gz"
    fastq_R2 = "${sample_name}_R2.fastq.gz"

    """
    # Convert to pairs of fastq files
    seqtk seq -1 ${fastq} | bgzip -@ ${task.cpus} > ${fastq_R1}
    seqtk seq -1 ${fastq} | bgzip -@ ${task.cpus} > ${fastq_R2}

    # Preprocess fasta file
    gzip -cd --fast ${fasta} > ${fasta_raw}
    bgzip -@ ${task.cpus} -c ${fasta_raw} > ${fasta_bgzip}
    samtools faidx ${fasta_bgzip}
    """
    
}