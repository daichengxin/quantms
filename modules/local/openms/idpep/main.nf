// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IDPEP {
    label 'process_very_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "openms::openms=2.7.0pre" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0"
    } else {
        container "quay.io/biocontainers/openms:2.6.0--h4afb90d_0"
    }

    input:
    tuple val(meta), path(id_file)

    output:
    tuple val(meta), path("${id_file.baseName}_fdr.idXML"), val("q-value_score"), emit: id_files_ForIDPEP
    path "*.version.txt", emit: version
    path "*.log", emit: log

    script:
    def software = getSoftwareName(task.process)

    """
    IDPosteriorErrorProbability \\
        -in ${id_file} \\
        -out ${id_file.baseName}_idpep.idXML \\
        -fit_algorithm:outlier_handling $params.outlier_handling \\
        -threads ${task.cpus} \\
        -debug 100 \\
        $options.args \\
        > ${id_file.baseName}_idpep.log

    echo \$(IDPosteriorErrorProbability 2>&1) > ${software}.version.txt
    """
}
