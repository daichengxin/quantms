
include { MSGF_DB_INDEXING } from '../../../modules/local/utils/msgf_db_indexing/main'
include { MSGF  } from '../../../modules/local/openms/msgf/main'
include { COMET } from '../../../modules/local/openms/comet/main'
include { SAGE  } from '../../../modules/local/openms/sage/main'
include { PSM_CLEAN            } from '../../../modules/local/utils/psm_clean/main'
include { MSRESCORE_FINE_TUNING} from '../../../modules/local/utils/msrescore_fine_tuning/main'
include { MSRESCORE_FEATURES   } from '../../../modules/local/utils/msrescore_features/main'

workflow PEPTIDE_DATABASE_SEARCH {
    take:
    ch_mzmls_search
    ch_searchengine_in_db

    main:
    (ch_id_msgf, ch_id_comet, ch_id_sage, ch_versions) = [ Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty() ]

    if (params.search_engines.contains("msgf")) {
        MSGF_DB_INDEXING(ch_searchengine_in_db)
        ch_versions = ch_versions.mix(MSGF_DB_INDEXING.out.versions)

        MSGF(ch_mzmls_search.combine(ch_searchengine_in_db).combine(MSGF_DB_INDEXING.out.msgfdb_idx))
        ch_versions = ch_versions.mix(MSGF.out.versions)
        ch_id_msgf = ch_id_msgf.mix(MSGF.out.id_files_msgf)
    }

    if (params.search_engines.contains("comet")) {
        COMET(ch_mzmls_search.combine(ch_searchengine_in_db))
        ch_versions = ch_versions.mix(COMET.out.versions)
        ch_id_comet = ch_id_comet.mix(COMET.out.id_files_comet)
    }

    // sorted mzmls to generate same batch ids when enable cache
    ch_mzmls_sorted_search = ch_mzmls_search.collect(flat: false, sort: { a, b -> a[0]["mzml_id"] <=> b[0]["mzml_id"] }).flatMap()
    if (params.search_engines.contains("sage")) {
        cnt = 0
        ch_meta_mzml_db = ch_mzmls_sorted_search.map{ metapart, mzml ->
            cnt++
            def groupkey = metapart.labelling_type +
                    metapart.dissociationmethod +
                    metapart.fixedmodifications +
                    metapart.variablemodifications +
                    metapart.precursormasstolerance +
                    metapart.precursormasstoleranceunit +
                    metapart.fragmentmasstolerance +
                    metapart.fragmentmasstoleranceunit +
                    metapart.enzyme
            // TODO this only works if the metakeys are all the same
            //  otherwise we need to group by key first and then batch
            def batch = cnt % params.sage_processes
            // TODO hash the key to make it shorter?
            [groupkey, batch, metapart, mzml]
        }
        // group into chunks to be processed at the same time on the same node by sage
        // TODO I guess if we parametrize the nr of files per process, it is more
        //  efficient (because this process can start as soon as this number of files
        //  are available and does not need to wait and see how many Channel entries
        //  belong to batch X). But the problem is groupTuple(size:) cannot be
        //  specified with an output from a Channel. The only way would be to,
        //  IN THE VERY BEGINNING, parse
        //  the number of files (=lines?) in the SDRF/design (outside of a process),
        //  save this value and pass it along the pipeline.
        ch_meta_mzml_db_chunked = ch_meta_mzml_db.groupTuple(by: [0,1])

        SAGE(ch_meta_mzml_db_chunked.combine(ch_searchengine_in_db))
        ch_versions = ch_versions.mix(SAGE.out.versions)
        // we can safely use merge here since it is the same process
        ch_id_sage = ch_id_sage.mix(SAGE.out.id_files_sage.transpose())
    }

    if (params.skip_rescoring != true && params.ms2features_enable == true) {
        // Only add ms2_model_dir if it's actually set and not empty
        // Handle cases where parameter might be empty string, null, boolean true, or whitespace
        // When --ms2features_model_dir is passed with no value, Nextflow may set it to boolean true
        if (params.ms2features_model_dir && params.ms2features_model_dir != true) {
            ms2_model_dir = Channel.from(file(params.ms2features_model_dir, checkIfExists: true))
        } else {
            ms2_model_dir = Channel.from(file("pretrained_models"))
        }

        if (params.ms2features_fine_tuning == true) {
            if (params.ms2features_generators.toLowerCase().contains('ms2pip')) {
                exit(1, 'Error: Fine tuning only supports AlphaPeptdeep!')
            } else {

                // Preparing train datasets and fine tuning MS2 model
                sage_train_datasets = ch_id_sage
                    .combine(ch_mzmls_search, by: 0)
                    .randomSample(params.fine_tuning_sample_run, 2025)
                    .combine(Channel.value("sage"))
                    .groupTuple(by: 3)

                msgf_train_datasets = ch_id_msgf
                    .combine(ch_mzmls_search, by: 0)
                    .randomSample(params.fine_tuning_sample_run, 2025)
                    .combine(Channel.value("msgf"))
                    .groupTuple(by: 3)

                comet_train_datasets = ch_id_comet
                    .combine(ch_mzmls_search, by: 0)
                    .randomSample(params.fine_tuning_sample_run, 2025)
                    .combine(Channel.value("comet"))
                    .groupTuple(by: 3)

                sage_train_datasets.mix(msgf_train_datasets)
                    .mix(comet_train_datasets)
                    .combine(ms2_model_dir)
                    .set { train_datasets }
                MSRESCORE_FINE_TUNING(train_datasets)
                ch_software_versions = ch_software_versions.mix(MSRESCORE_FINE_TUNING.out.versions)

                sage_features_input = Channel.value("sage").combine(ch_id_files_branched.sage)
                    .combine(MSRESCORE_FINE_TUNING.out.model_weight, by:0)
                msgf_features_input = Channel.value("msgf").combine(ch_id_files_branched.msgf)
                    .combine(MSRESCORE_FINE_TUNING.out.model_weight, by:0)
                comet_features_input = Channel.value("comet").combine(ch_id_files_branched.comet)
                    .combine(MSRESCORE_FINE_TUNING.out.model_weight, by:0)
                sage_features_input.mix(msgf_features_input).mix(comet_features_input)
                    .map { [it[1], it[2], it[3], it[4]] }
                    .set { ch_features_input }

                MSRESCORE_FEATURES(ch_features_input)
                ch_software_versions = ch_software_versions.mix(MSRESCORE_FEATURES.out.versions)
                ch_id_files_feats = MSRESCORE_FEATURES.out.idxml
            }
        } else{
            MSRESCORE_FEATURES(ch_id_files.combine(ch_file_preparation_results, by: 0).combine(ms2_model_dir))
            ch_software_versions = ch_software_versions.mix(MSRESCORE_FEATURES.out.versions)
            ch_id_files_feats = MSRESCORE_FEATURES.out.idxml
        }

    } else if (params.psm_clean == true) {
        ch_id_files = ch_id_msgf.mix(ch_id_comet).mix(ch_id_sage)
        PSM_CLEAN(ch_id_files.combine(ch_file_preparation_results, by: 0))
        ch_id_files_feats = PSM_CLEAN.out.idxml
        ch_software_versions = ch_software_versions.mix(PSM_CLEAN.out.versions)
    }

    emit:
    ch_id_files_idx = ch_id_files_feats
    versions        = ch_versions
}
