process ANALYZE_GENE_COVERAGE {
    tag "$meta.id"
    label 'process_single'
    label 'publish_outdir'

    input:
    tuple val(meta), path(deduped_bam), path(deduped_bai)
    path(genes_bed)

    output:
    tuple val(meta), path(genes_coverage), emit: genes_coverage

    script:
    def id          = meta.id
    def fasta       = meta.fasta

    genes_coverage  = "${id}.coverage.txt"

    """
    bedtools coverage -a ${genes_bed} -b ${deduped_bam} > ${genes_coverage}
    """
}