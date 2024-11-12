process VARIANT_CALLING {
    tag "$meta.id"
    label 'process_medium'
    label 'publish_outdir'

    input:
    tuple val(meta), path(deduped_bam), path(deduped_bai)
    path(fasta_dir)

    output:
    tuple val(meta), path(vcf), path(vcf_gz), path(vcf_tbi), emit: vcf

    script:
    def id          = meta.id
    def fasta       = meta.fasta
    def vcf_hc      = "hc.vcf.gz"

    vcf             = "${id}.vcf"
    vcf_gz          = "${vcf}.gz"
    vcf_tbi         = "${vcf_gz}.tbi"

    """
    # Variant Haplotype Caller
    gatk --java-options "-Xmx16g" HaplotypeCallerSpark  \
    -R ${fasta} \
    -I ${deduped_bam} \
    -O ${vcf_hc} \
    --sample-name ${id} \
    -ERC GVCF

    # Perform tabix & bgzip
    zless ${vcf_hc} > ${vcf}
    bgzip -cf ${vcf} > ${vcf_gz}
    tabix ${vcf_gz}
    """
}