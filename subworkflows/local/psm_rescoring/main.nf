//
// Extract psm feature and ReScoring psm
//

include { PERCOLATOR              } from '../../../modules/local/openms/percolator/main'
include { MSRESCORE_FEATURES      } from '../../../modules/local/utils/msrescore_features/main'
include { GET_SAMPLE              } from '../../../modules/local/utils/extract_sample/main'
include { ID_MERGER               } from '../../../modules/local/openms/id_merger/main'
include { ID_RIPPER               } from '../../../modules/local/openms/id_ripper/main'
include { SPECTRUM_FEATURES       } from '../../../modules/local/utils/spectrum_features/main'
include { PSM_CLEAN               } from '../../../modules/local/utils/psm_clean/main'
include { MSRESCORE_FINE_TUNING   } from '../../../modules/local/utils/msrescore_fine_tuning/main'

workflow PSM_RESCORING {
    take:
    ch_file_preparation_results
    ch_id_files
    ch_expdesign

    main:
    ch_software_versions = Channel.empty()
    ch_results  = Channel.empty()
    ch_fdridpep = Channel.empty()

    if (params.ms2features_enable == true) {

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
                // Split ch_id_files by search_engines
                ch_id_files.combine(ch_file_preparation_results, by: 0).branch{ meta, filename, mzml_name ->
                    sage: filename.name.contains('sage')
                        return [meta, filename, mzml_name]
                    msgf: filename.name.contains('msgf')
                        return [meta, filename, mzml_name]
                    comet: filename.name.contains('comet')
                        return [meta, filename, mzml_name]
                }.set{ch_id_files_branched}

                // Preparing train datasets and fine tuning MS2 model
                sage_train_datasets = ch_id_files_branched.sage.randomSample(params.fine_tuning_sample_run, 2025).combine(
                    Channel.value("sage")
                ).groupTuple(by: 3)
                msgf_train_datasets = ch_id_files_branched.msgf.randomSample(params.fine_tuning_sample_run, 2025).combine(
                    Channel.value("msgf")
                ).groupTuple(by: 3)
                comet_train_datasets = ch_id_files_branched.comet.randomSample(params.fine_tuning_sample_run, 2025).combine(
                    Channel.value("comet")
                ).groupTuple(by: 3)
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
        } else{
            MSRESCORE_FEATURES(ch_id_files.combine(ch_file_preparation_results, by: 0).combine(ms2_model_dir))
            ch_software_versions = ch_software_versions.mix(MSRESCORE_FEATURES.out.versions)
            ch_id_files_feats = MSRESCORE_FEATURES.out.idxml
        }

    } else if (params.psm_clean == true) {
        PSM_CLEAN(ch_id_files.combine(ch_file_preparation_results, by: 0))
        ch_id_files_feats = PSM_CLEAN.out.idxml
        ch_software_versions = ch_software_versions.mix(PSM_CLEAN.out.versions)
    } else {
        ch_id_files_feats = ch_id_files
    }

    // Add SNR features to percolator
    if (params.ms2features_snr) {
        SPECTRUM_FEATURES(ch_id_files_feats.combine(ch_file_preparation_results, by: 0))
        ch_id_files_feats = SPECTRUM_FEATURES.out.id_files_snr
        ch_software_versions = ch_software_versions.mix(SPECTRUM_FEATURES.out.versions)
    }

    // Rescoring for independent run, Sample or whole experiments
    if (params.ms2features_range == "independent_run") {
        PERCOLATOR(ch_id_files_feats)
        ch_software_versions = ch_software_versions.mix(PERCOLATOR.out.versions)
        ch_consensus_input = PERCOLATOR.out.id_files_perc
    } else if (params.ms2features_range == "by_sample") {
        // Sample map
        GET_SAMPLE(ch_expdesign)
        ch_software_versions = ch_software_versions.mix(GET_SAMPLE.out.versions)
        ch_expdesign_sample = GET_SAMPLE.out.ch_expdesign_sample
        ch_expdesign_sample.splitCsv(header: true, sep: '\t')
            .map { get_sample_map(it) }.set{ sample_map_idv }

        ch_id_files_feats.map {[it[0].mzml_id, it[0], it[1]]}.set { ch_id_files_feats}
        ch_id_files_feats.combine(sample_map_idv, by: 0).map {[it[1], it[2], it[3]]}.set{ch_id_files_feats}

        // Group by search_engines and convert meta
        ch_id_files_feats.branch{ meta, filename, sample  ->
            sage: filename.name.contains('sage')
                return [meta, filename, sample]
            msgf: filename.name.contains('msgf')
                return [meta, filename, sample]
            comet: filename.name.contains('comet')
                return [meta, filename, sample]
        }.set{ch_id_files_feat_branched}

        // ID_MERGER for samples group
        ID_MERGER(ch_id_files_feat_branched.comet.groupTuple(by: 2)
            .mix(ch_id_files_feat_branched.msgf.groupTuple(by: 2))
            .mix(ch_id_files_feat_branched.sage.groupTuple(by: 2)))
        ch_software_versions = ch_software_versions.mix(ID_MERGER.out.versions)

        PERCOLATOR(ID_MERGER.out.id_merged)
        ch_software_versions = ch_software_versions.mix(PERCOLATOR.out.versions)

        // Currently only ID runs on exactly one mzML file are supported in CONSENSUSID. Split idXML by runs
        ID_RIPPER(PERCOLATOR.out.id_files_perc)
        ch_file_preparation_results.map{[it[0].mzml_id, it[0]]}.set{meta}
        ID_RIPPER.out.id_rippers.flatten().map { add_file_prefix (it)}.set{id_rippers}
        meta.combine(id_rippers, by: 0)
                .map{ [it[1], it[2]]}
                .set{ ch_consensus_input }
        ch_software_versions = ch_software_versions.mix(ID_RIPPER.out.versions)

    } else if (params.ms2features_range == "by_project"){
        ch_id_files_feats.map {[it[0].experiment_id, it[0], it[1]]}.set { ch_id_files_feats}

        // Split ch_id_files_feats by search_engines
        ch_id_files_feats.branch{ experiment_id, meta, filename ->
            sage: filename.name.contains('sage')
                return [meta, filename, experiment_id]
            msgf: filename.name.contains('msgf')
                return [meta, filename, experiment_id]
            comet: filename.name.contains('comet')
                return [meta, filename, experiment_id]
        }.set{ch_id_files_feat_branched}

        // ID_MERGER for whole experiments
        ID_MERGER(ch_id_files_feat_branched.comet.groupTuple(by: 2)
            .mix(ch_id_files_feat_branched.msgf.groupTuple(by: 2))
            .mix(ch_id_files_feat_branched.sage.groupTuple(by: 2)))
        ch_software_versions = ch_software_versions.mix(ID_MERGER.out.versions)

        PERCOLATOR(ID_MERGER.out.id_merged)
        ch_software_versions = ch_software_versions.mix(PERCOLATOR.out.versions)

        // Currently only ID runs on exactly one mzML file are supported in CONSENSUSID. Split idXML by runs
        ID_RIPPER(PERCOLATOR.out.id_files_perc)
        ch_file_preparation_results.map{[it[0].mzml_id, it[0]]}.set{meta}
        ID_RIPPER.out.id_rippers.flatten().map { add_file_prefix (it)}.set{id_rippers}
        meta.combine(id_rippers, by: 0)
                .map{ [it[1], it[2]]}
                .set{ ch_consensus_input }
        ch_software_versions = ch_software_versions.mix(ID_RIPPER.out.versions)
    }
    ch_rescoring_results = ch_consensus_input

    emit:
    results = ch_rescoring_results
    versions = ch_software_versions
}

def add_file_prefix(file_path) {
    position = file(file_path).name.lastIndexOf('_sage_perc.idXML')
    if (position == -1) {
        position = file(file_path).name.lastIndexOf('_comet_perc.idXML')
        if (position == -1) {
            position = file(file_path).name.lastIndexOf('_msgf_perc.idXML')
        }
    }
    file_name = file(file_name).name.take(position)
    return [file_name, file_path]
}

// Function to get sample map
def get_sample_map(LinkedHashMap row) {
    def sample_map = [:]

    filestr               = row.Spectra_Filepath
    file_name             = file(filestr).name.take(file(filestr).name.lastIndexOf('.'))
    sample                = row.Sample

    return [file_name, sample]

}
