process PERCOLATOR {
    tag "$meta.mzml_id"
    label 'process_medium'
    label 'openms'

    conda "bioconda::openms-thirdparty=3.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:3.2.0--h9ee0642_4' :
        'biocontainers/openms-thirdparty:3.2.0--h9ee0642_4' }"

    input:
    tuple val(meta), path(id_file)

    output:
    tuple val(meta), path("*_perc.idXML"), val("MS:1001491"), emit: id_files_perc
    path "versions.yml", emit: versions
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"

    """
    OMP_NUM_THREADS=$task.cpus PercolatorAdapter \\
        -in ${id_file} \\
        -out ${id_file.baseName}_perc.idXML \\
        -threads $task.cpus \\
        -subset_max_train $params.subset_max_train \\
        -decoy_pattern $params.decoy_string \\
        -post_processing_tdc \\
        -score_type pep \\
        $args \\
        2>&1 | tee ${id_file.baseName}_percolator.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PercolatorAdapter: \$(PercolatorAdapter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
        percolator: \$(percolator -h 2>&1 | grep -E '^Percolator version(.*)' | sed 's/Percolator version //g')
    END_VERSIONS
    """
}
