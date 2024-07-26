process IDFILECONVERTER {
    label 'process_low'
    label 'openms'

    conda "bioconda::openms-thirdparty=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:3.1.0--h9ee0642_1' :
        'biocontainers/openms-thirdparty:3.1.0--h9ee0642_1' }"

    input:
    tuple val(meta), path(id_file), path(mz_file)

    output:
    tuple val(meta), path("${id_file.baseName}.pepXML"), emit: pepXML
    path "versions.yml", emit: version
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''

    """
    IDFileConverter \\
        -in ${id_file} \\
        -mz_file ${mz_file} \\
        -threads $task.cpus \\
        -out ${id_file.baseName}.pepXML \\
        $args \\
        2>&1 | tee ${id_file.baseName}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IDFileConverter: \$(IDFileConverter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
    END_VERSIONS
    """

}
