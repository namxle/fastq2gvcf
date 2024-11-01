/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
checkPathParamList = [
   params.fastq,
   params.fasta
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

ch_fasta = Channel.fromPath(params.fasta)
ch_fastq = Channel.fromPath(params.fastq)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { PREPROCESS        } from '../modules/preprocess'
include { ALIGNMENT_SORTING } from '../modules/alignment_sorting'


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


workflow ANALYZE {
    ch_input = create_input_channel(ch_fastq)

    PREPROCESS(
        ch_input,
        ch_fasta
    )

    ALIGNMENT_SORTING(
        PREPROCESS.out.fastq,
        PREPROCESS.out.fasta
    )
}

def create_input_channel (ch_fastq) {
    // create meta map
    def meta = [:]
    meta.id                 = params.sample_id
    meta.name               = params.sample_name

    return Channel.of([meta]).combine(ch_fastq)
}
