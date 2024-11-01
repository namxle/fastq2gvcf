process PREPROCESS {
    tag "$meta.id"
    label 'process_fixed_medium_cpu'
    label 'publish_outdir'

    input:
    tuple val(meta), path(fastq)
    path(fasta)

    output:
    tuple val(meta), path(fastq_R1), path(fastq_R2), emit: fastq
    tuple path(fasta), path(fasta_index)

    script:
    def sample_name = meta.name

    fasta_index = "${fasta}.fai"

    fastq_R1 = "${sample_name}_R1.fastq.gz"
    fastq_R2 = "${sample_name}_R2.fastq.gz"

    """
    seqtk seq -1 ${fastq} | bgzip -@ ${task.cpus} > ${fastq_R1}
    seqtk seq -1 ${fastq} | bgzip -@ ${task.cpus} > ${fastq_R2}

    samtools faidx ${fasta}
    """
    
}