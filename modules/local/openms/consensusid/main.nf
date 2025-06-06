process CONSENSUSID {
    tag "$meta.mzml_id"
    label 'process_single'
    label 'openms'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/bigbio/openms-tools-thirdparty-sif:2025.04.14' :
        'ghcr.io/bigbio/openms-tools-thirdparty:2025.04.14' }"

    input:
    tuple val(meta), path(id_file), val(qval_score)

    output:
    tuple val(meta), path("${meta.mzml_id}_consensus.idXML"), emit: consensusids
    path "versions.yml", emit: versions
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"

    """
    ConsensusID \\
        -in ${id_file} \\
        -out ${meta.mzml_id}_consensus.idXML \\
        -per_spectrum \\
        -threads $task.cpus \\
        -algorithm $params.consensusid_algorithm \\
        -filter:min_support $params.min_consensus_support \\
        -filter:considered_hits $params.consensusid_considered_top_hits \\
        -debug $params.consensusid_debug \\
        $args \\
        2>&1 | tee ${meta.mzml_id}_consensusID.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ConsensusID: \$(ConsensusID 2>&1  | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
    END_VERSIONS
    """
}
