process MARK_DUPLICATES {
    tag "$meta.id"
    label 'process_dual'
    label 'publish_outdir'

    input:
    tuple val(meta), path(bam), path(bam_bai)
    path(fasta_dir)

    output:
    tuple val(meta), path(deduped_bam), path(deduped_bai), emit: deduped_bam

    script:
    def id          = meta.id
    def fasta       = meta.fasta

    deduped_bam     = "${id}.deduped.bam"
    deduped_bai     = "${deduped_bam}.bai"


    """
    # Mark duplicate 
    gatk MarkDuplicatesSpark \
      -I ${bam} \
      -O ${deduped_bam} \
      --remove-all-duplicates true
      
    # Index deduped bam
    samtools index ${deduped_bam}
    """
}