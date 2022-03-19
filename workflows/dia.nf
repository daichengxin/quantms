/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES: Local to the pipeline
//
include { DIANNSEARCH } from '../modules/local/diannsearch/main'
include { GENERATE_DIANN_CFG  as DIANNCFG } from '../modules/local/generate_diann_cfg/main'
include { CONVERT2MSSTATS } from '../modules/local/convert2msstats/main'
include { LIBRARYGENERATION } from '../modules/local/librarygeneration/main'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DIA {
    take:
    file_preparation_results
    ch_expdesign

    main:

    ch_software_versions = Channel.empty()

    file_preparation_results.multiMap {
                                meta: it[0]
                                mzml: it[1]
                                }
                            .set { result }

    DIANNCFG(result.meta.collect(), result.mzml.collect())
    ch_software_versions = ch_software_versions.mix(DIANNCFG.out.version.ifEmpty(null))

    LIBRARYGENERATION(Channel.fromPath(params.database), DIANNCFG.out.library_config)

    DIANNSEARCH(DIANNCFG.out.mzmls_for_diann.collect(), LIBRARYGENERATION.out.lib_splib, DIANNCFG.out.search_cfg)
    ch_software_versions = ch_software_versions.mix(DIANNSEARCH.out.version.ifEmpty(null))

    CONVERT2MSSTATS(DIANNSEARCH.out.report, ch_expdesign)
    versions        = ch_software_versions

    emit:
    versions    = versions
}

/*
========================================================================================
    THE END
========================================================================================
*/
