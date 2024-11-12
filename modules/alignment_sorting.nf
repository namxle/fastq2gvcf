process ALIGNMENT_SORTING {
    tag "$meta.id"
    label 'process_fixed_medium_cpu'
    label 'publish_outdir'

    input:
    tuple val(meta), path(fastq_R1), path(fastq_R2)
    path(fasta_data)

    output:
    tuple val(meta), path(sorted_bam), path(sorted_bam_bai), emit: bam

    script:
    def bwa_bam     = "bwa.bam"
    def ali_bam     = "ali.bam"
    def fasta       = meta.fasta
    def delimiter   = "\\t"
    def platform    = "ILLUMINA"    
    def group       = "${meta.name}"
    def sample      = "${meta.name}"

    sorted_bam      = "sorted.bam"
    sorted_bam_bai  = "${sorted_bam}.bai"

    """
    # Alignment
    bwa mem \
    -t ${task.cpus} \
    -M \
    -R \
    "@RG${delimiter}ID:${group}${delimiter}SM:${sample}${delimiter}PL:${platform}${delimiter}LB:default" \
    ${fasta} \
    ${fastq_R1} \
    ${fastq_R2} | \
    samtools view -Sb -o ${bwa_bam} -@ ${task.cpus}

    # Remove un-mapped reads
    samtools view -F 4 ${bwa_bam} -o ${ali_bam} --threads ${task.cpus}
    samtools sort ${ali_bam} -o ${sorted_bam}
    samtools index ${sorted_bam}
    """
    
}