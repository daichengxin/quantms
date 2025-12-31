process DOWNLOAD_MODEL {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/daichengxin/quantms-rescoring-sif:0.0.14' :
        'ghcr.io/daichengxin/quantms-rescoring:0.0.14' }"

    input:
    val(model_list)

    output:
    path "rescore_model" , emit: model_weights
    path "versions.yml"     , emit: versions
    path "*.log"            , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    rescoring download_models \\
        --models ${model_list} \\
        --model_dir ./rescore_model \\
        $args \\
        2>&1 | tee download_models.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quantms-rescoring: \$(rescoring --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+')
        ms2pip: \$(ms2pip --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+')
        deeplc: \$(deeplc --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+')
        MS2Rescore: \$(ms2rescore --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n 1)
    END_VERSIONS
    """
}
