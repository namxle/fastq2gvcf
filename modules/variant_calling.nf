process VARIANT_CALLING {
    tag "$meta.id"
    label 'process_fixed_medium_cpu'
    label 'publish_outdir'

    input:
    tuple val(meta), path(bam)
    path(fasta_dir)

    output:
    tuple val(meta), path(vcf), emit: vcf

    script:
    def fasta   = meta.fasta
    def id      = meta.id

    vcf         = "${id}.vcf.gz"

    """
    gatk --java-options "-Xmx4g" HaplotypeCaller  \
    -R ${fasta} \
    -I ${bam} \
    -O ${vcf} \
    -ERC GVCF
    """
}