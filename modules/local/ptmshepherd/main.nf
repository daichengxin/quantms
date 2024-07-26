process PTMSHEPHERD {
    tag "$meta.mzml_id"
    label 'process_medium'

    conda "bioconda::fragpipe=20.0"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fragpipe:20.0--hdfd78af_3"
    } else {
        container "biocontainers/fragpipe:20.0--hdfd78af_3"
    }

    input:
    tuple val(meta), path(psm_tsv), path(mzML), path(ptmshepherd_cli), path(fasta), path(config)

    output:
    path "global.modsummary.tsv", emit: psm_info
    path "global.profile.tsv", emit: pepXML
    path "versions.yml", emit: version
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"


    """
    java -jar "${ptmshepherd-2.0.5_CLI.jar}" \\
        ${config} \\
        2>&1 | tee ptmshepherd.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PTMShepherd: 2.0.5')
    END_VERSIONS
    """

}
