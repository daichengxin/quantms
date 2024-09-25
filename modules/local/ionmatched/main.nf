process IONMATCHED {
    tag "$csv_file.baseName"
    label 'process_medium'

    conda "bioconda::spectrum_utils=0.4.2"
    container "docker.io/daicx1/spectrum_annotations:v1.0"
    containerOptions = "--user root"

    input:
    tuple val(meta), path(csv_file)

    output:
    path "*_psm.parquet", emit: psm_info
    path "versions.yml", emit: version
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''


    """
    ions_annotation.py "${csv_file}" \\
        $meta.fragmentmasstolerance \\
        $meta.fragmentmasstoleranceunit \\
        2>&1 | tee ions_annotation.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyopenms: \$(pip show pyopenms | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """
}
