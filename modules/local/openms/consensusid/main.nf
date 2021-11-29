// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CONSENSUSID {
    label 'process_medium'
    // TODO could be easily parallelized
    label 'process_single_thread'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::gnuplot openms::openms-thirdparty=2.7.0pre" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://ftp.pride.ebi.ac.uk/pride/data/tools/quantms-dev.sif"
    } else {
        container "quay.io/bigbio/quantms:dev"
    }

    input:
    tuple val(meta), path(id_file), val(qval_score)

    output:
    tuple val(meta), path("${meta.id}_consensus.idXML"), emit: consensusids
    path "*.version.txt", emit: version
    path "*.log", emit: log
    script:
    def software = getSoftwareName(task.process)

    """
    ConsensusID \\
        -in ${id_file} \\
        -out ${meta.id}_consensus.idXML \\
        -per_spectrum \\
        -threads $task.cpus \\
        -algorithm $params.consensusid_algorithm \\
        -filter:min_support $params.min_consensus_support \\
        -filter:considered_hits $params.consensusid_considered_top_hits \\
        -debug 100 \\
        $options.args \\
        > ${meta.id}_consensusID.log

    echo \$(ConsensusID 2>&1) > ${software}.version.txt
    """
}
