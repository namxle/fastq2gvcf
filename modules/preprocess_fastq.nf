process PREPROCESS_FASTQ {
    tag "$meta.id"
    label 'process_medium'
    label 'publish_outdir'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path(fastq_R1), path(fastq_R2), emit: fastq

    script:
    def sample_name = meta.name
    def fasta       = "."
    def fasta_raw   = "."

    fastq_R1 = "${sample_name}_R1.fastq.gz"
    fastq_R2 = "${sample_name}_R2.fastq.gz"

    """
    # Convert to pairs of fastq files
    seqtk seq -1 ${fastq} | bgzip -@ ${task.cpus} > ${fastq_R1}
    seqtk seq -2 ${fastq} | bgzip -@ ${task.cpus} > ${fastq_R2}

    # Preprocess fasta file
    # gzip -cd --fast ${fasta} > ${fasta_raw}
    # bwa index ${fasta_raw}
    """
    
}