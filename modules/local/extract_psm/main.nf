process PSMCONVERSION {
    tag "$meta.mzml_id"
    label 'process_medium'

    conda "bioconda::pmultiqc=0.0.25"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pmultiqc:0.0.25--pyhdfd78af_0' :
        'biocontainers/pmultiqc:0.0.25--pyhdfd78af_0' }"


    input:
    tuple val(meta), path(idxml_file), path(spectrum_df), path(exp_design)

    output:
    path "*_psm.csv", emit: psm_info
    path "versions.yml", emit: version
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"


    """
    psm_conversion.py "${idxml_file}" \\
        ${spectrum_df} \\
        ${exp_design} \\
        $params.enable_mod_localization \\
        $params.export_decoy_psm \\
        2>&1 | tee extract_idxml.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyopenms: \$(pip show pyopenms | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """
}
