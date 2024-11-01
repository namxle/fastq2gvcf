process ALIGNMENT_SORTING {
    tag "$meta.id"
    label 'process_fixed_medium_cpu'
    label 'publish_outdir'

    input:
    tuple val(meta), path(fastq_R1), path(fastq_R2)
    tuple path(fasta), path(fasta_index)

    output:
    tuple val(meta), path(sorted_bam)

    script:
    def sample_name = meta.name
    def bwa_bam     = "bwa.bam"
    def ali_bam     = "ali.bam"

    sorted_bam      = "sorted.bam"

    """
    # Align
    bwa mem \
    -t ${task.cpus} \
    -M \
    ${fasta} \
    ${fastq_R1} \
    ${fastq_R2} | \
    samtools view -Sb -o ${bwa_bam} -@ ${task.cpus}

    # Remove un-mapped reads
    samtools view -F 4 ${bwa_bam} -o ${ali_bam} --threads ${task.cpus}
    samtools sort ${ali_bam} -o ${sorted_bam}
    """
    
}