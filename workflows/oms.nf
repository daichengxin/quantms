//
// MODULE: Local to the pipeline
//
include { DECOYDATABASE } from '../modules/local/openms/decoydatabase/main'
include { CONSENSUSID   } from '../modules/local/openms/consensusid/main'
include { EXPORTPSMTSV  } from '../modules/local/openms/exportpsmtsv/main'
include { IDFILECONVERTER } from '../modules/local/openms/idfileconverter/main'
include { IDMERGER           } from '../modules/local/openms/idmerger/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { DATABASESEARCHENGINES } from '../subworkflows/local/databasesearchengines'
include { PSMRESCORING          } from '../subworkflows/local/psmrescoring'
include { PSMFDRCONTROL         } from '../subworkflows/local/psmfdrcontrol'
include { PHOSPHOSCORING        } from '../subworkflows/local/phosphoscoring'
include { PROTEININFERENCE } from '../subworkflows/local/proteininference'


workflow OMS {
    take:
    ch_file_preparation_results
    ch_database_wdecoy
    ch_spectrum_data
    ch_expdesign

    main:

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: DatabaseSearchEngines
    //
    DATABASESEARCHENGINES (
        ch_file_preparation_results,
        ch_database_wdecoy
    )
    ch_software_versions = ch_software_versions.mix(DATABASESEARCHENGINES.out.versions.ifEmpty(null))

    //
    // SUBWORKFLOW: PSMReScoring
    //
    PSMRESCORING (ch_file_preparation_results, DATABASESEARCHENGINES.out.ch_id_files_idx, ch_expdesign)
    ch_software_versions = ch_software_versions.mix(PSMRESCORING.out.versions.ifEmpty(null))

    //
    // SUBWORKFLOW: PSMFDRCONTROL
    //
    ch_psmfdrcontrol     = Channel.empty()
    ch_consensus_results = Channel.empty()
    if (params.search_engines.split(",").size() > 1) {
        CONSENSUSID(PSMRESCORING.out.results.groupTuple(size: params.search_engines.split(",").size()))
        ch_software_versions = ch_software_versions.mix(CONSENSUSID.out.version.ifEmpty(null))
        ch_psmfdrcontrol = CONSENSUSID.out.consensusids
        ch_consensus_results = CONSENSUSID.out.consensusids
    } else {
        ch_psmfdrcontrol = PSMRESCORING.out.results
    }

    PSMFDRCONTROL(ch_psmfdrcontrol)
    ch_software_versions = ch_software_versions.mix(PSMFDRCONTROL.out.version.ifEmpty(null))

    //
    // SUBWORKFLOW：PHOSPHOSCORING
    //
    if (params.enable_mod_localization) {
        PHOSPHOSCORING(ch_file_preparation_results, PSMFDRCONTROL.out.id_filtered)
        ch_software_versions = ch_software_versions.mix(PHOSPHOSCORING.out.version.ifEmpty(null))
        ch_id_results = PHOSPHOSCORING.out.id_luciphor
    } else {
        ch_id_results = PSMFDRCONTROL.out.id_filtered
    }

    // 蛋白质推断打分：两种选择 ProteinInference, TPP ProteinProphet
    //
    // MODULE: FILEMERGE
    //
    IDMERGER(PSMFDRCONTROL.out.id_filtered.collect())
    ch_software_versions = ch_software_versions.mix(IDMERGER.out.version.ifEmpty(null))

    //
    // SUBWORKFLOW: PROTEININFERENCE
    //
    PROTEININFERENCE(IDMERGER.out.id_merge)
    ch_software_versions = ch_software_versions.mix(PROTEININFERENCE.out.version.ifEmpty(null))

    // idXML2pepXML    iProphet和Percolator等价
    IDFILECONVERTER(PSMFDRCONTROL.out.id_filtered.combine(ch_file_preparation_results, by: 0))

    // exportPSMTSV
    EXPORTPSMTSV(IDFILECONVERTER.out.pepXML.combine(PSMFDRCONTROL.out.id_filtered, by: 0))

    // PTMProphet  修饰定位打分
    PTMPROPHET(EXPORTPSMTSV.out.pepXML)


    emit:
    id_results              = ch_id_results
    psmrescoring_results    = PSMRESCORING.out.results
    ch_consensus_results    = ch_consensus_results
    version                 = ch_software_versions
}
