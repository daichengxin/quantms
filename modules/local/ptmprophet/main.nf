process PTMPROPHET {
    tag "$pepxml_file.baseName"
    label 'process_medium'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fragpipe:20.0--hdfd78af_3"
    } else {
        container "fcyucn/fragpipe:22.0"
    }

    input:
    path(pepxml_file)

    output:
    path("${pepxml_file.baseName}_ptmprophet.pepXML"), emit: ptmprophet_pepXML
    path "versions.yml", emit: version
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${pepxml_file.baseName}"


    """
    /fragpipe_bin/fragPipe-22.0/fragpipe/tools/PTMProphet
        "${pepxml_file}" \\
        EM=2 \\
        VERBOSE=true \\
        FRAGPPMTOL=$params.fragptmtol \\
        2>&1 | tee PTMProphet.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PTMProphet: \$(/fragpipe_bin/fragPipe-22.0/fragpipe/tools/PTMProphet | grep "VERSION" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """


}
