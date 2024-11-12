process VEP_ANNOTATION {
    tag "$meta.id"
    label 'process_fixed_medium_cpu'
    label 'publish_outdir'

    input:
    tuple val(meta), path(vcf), path(vcf_gz), path(vcf_tbi)

    output:
    tuple val(meta), path(vcf), emit: vcf

    script:
    def id          = meta.id
    def fasta       = meta.fasta
    def deduped_bam = "${id}.deduped.bam"

    vcf             = "${id}.vcf.gz"
}