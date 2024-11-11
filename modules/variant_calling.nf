process VARIANT_CALLING {
    tag "$meta.id"
    label 'process_fixed_medium_cpu'
    label 'publish_outdir'

    input:
    tuple val(meta), path(bam), path(bam_bai)
    path(fasta_dir)

    output:
    tuple val(meta), path(vcf), emit: vcf

    script:
    def id          = meta.id
    def fasta       = meta.fasta
    def deduped_bam = "${id}.deduped.bam"

    vcf             = "${id}.vcf.gz"

    """
    # Mark duplicate 
    gatk MarkDuplicatesSpark \
      -I ${bam} \
      -O ${deduped_bam} \
      --remove-all-duplicates true
      
    # Index deduped bam
    samtools index ${deduped_bam}

    # Variant Haplotype Caller
    gatk --java-options "-Xmx16g" HaplotypeCaller  \
    -R ${fasta} \
    -I ${deduped_bam} \
    -O ${vcf} \
    --sample-name ${id} \
    -ERC GVCF
    """
}