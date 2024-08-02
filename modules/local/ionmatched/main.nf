process IONMATCHED {
    tag "$csv_file.baseName"
    label 'process_medium'

    conda "bioconda::spectrum_utils=0.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spectrum_utils:0.4.2--pyhdfd78af_0' :
        'biocontainers/spectrum_utils:0.4.2--pyhdfd78af_0' }"


    input:
    path(csv_file)

    output:
    path "*_psm.csv", emit: psm_info
    path "versions.yml", emit: version
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''


    """
    ions_annotation.py "${csv_file}" \\
        2>&1 | tee ions_annotation.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyopenms: \$(pip show pyopenms | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """
}
