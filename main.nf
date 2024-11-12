/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fastq2gvcf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/namxle/fastq2gvcf
----------------------------------------------------------------------------------------
*/


nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

Workflow.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

def fasta_dir = "${projectDir}/assets/fasta"
def genes_bed = "${projectDir}/assets/genes.bed"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
checkPathParamList = [
   params.fastq,
   genes_bed
]

for (param in checkPathParamList) { 
    if (param) { 
        file_path = file(param)
        if (!file_path.exists()) { 
            println "File check: ${param} not exist!"
            NfcoreTemplate.workflowError(params)
        }
    } 
}

ch_fastq = Channel.fromPath(params.fastq)
ch_fasta = Channel.fromPath("${fasta_dir}/*", type: 'any', hidden: true)

ch_genes_bed = Channel.fromPath(genes_bed)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { PREPROCESS_FASTQ      } from './modules/preprocess_fastq'
include { ALIGNMENT_SORTING     } from './modules/alignment_sorting'
include { MARK_DUPLICATES       } from './modules/mark_duplicates'
include { ANALYZE_GENE_COVERAGE } from './modules/analyze_gene_coverage'

// include { VARIANT_CALLING   } from './modules/variant_calling'


//
// MODULE: Loaded from subworkflows/local-fastq/
//


//  
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow {
    ch_input = create_input_channel(ch_fastq)

    PREPROCESS_FASTQ(
        ch_input
    )

    ALIGNMENT_SORTING(
        PREPROCESS.out.fastq,
        ch_fasta.collect()
    )

    // VARIANT_CALLING(
    //     ALIGNMENT_SORTING.out.bam,
    //     ch_fasta.collect()
    // )

    MARK_DUPLICATES(
        ALIGNMENT_SORTING.out.bam,
        ch_fasta.collect()
    )

    ANALYZE_GENE_COVERAGE(
        MARK_DUPLICATES.out.deduped_bam,
        ch_genes_bed
    )
}


def create_input_channel (ch_fastq) {
    // create meta map
    def meta = [:]
    meta.id                 = params.sample_id
    meta.name               = params.sample_name
    meta.fasta              = "ref.fa"

    return Channel.of([meta]).combine(ch_fastq)
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    println "Execution status: ${ workflow.success ? 'OK' : 'FAILED' }"

    // Call api for success
    if (workflow.success) {
        // NfcoreTemplate.submit(params.submit_url, params.sample_id, 'success')
        println "Success"
    }
    println "Done"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onError {
    NfcoreTemplate.submit(params.submit_url, params.sample_id, 'error')
    println "Error"
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/