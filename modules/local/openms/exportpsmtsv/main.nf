process EXPORTPSMTSV {
    tag "$meta.mzml_id"
    label 'process_medium'

    conda "bioconda::pyopenms=3.1.0"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pyopenms:3.1.0--py39h9b8898c_0"
    } else {
        container "biocontainers/pyopenms:3.1.0--py39h9b8898c_0"
    }

    input:
    tuple val(meta), path(pepxml_file), path(idxml_file)

    output:
    path "psm.tsv", emit: psm_info
    path "corrected.pepXML", emit: pepXML
    path "versions.yml", emit: version
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"


    """
    exportpsm2tsv "${pepxml_file}" \\
        ${idxml_file} \\
        2>&1 | tee extract_idxml.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyopenms: \$(pip show pyopenms | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """

}
