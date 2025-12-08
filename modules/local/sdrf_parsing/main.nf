process SDRF_PARSING {
    tag "$sdrf.Name"
    label 'process_tiny'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quantms-utils:0.0.23--pyh7e72e81_0' :
        'biocontainers/quantms-utils:0.0.23--pyh7e72e81_0' }"

    input:
    path sdrf

    output:
    path "${sdrf.baseName}_openms_design.tsv", emit: ch_expdesign
    path "${sdrf.baseName}_config.tsv"       , emit: ch_sdrf_config_file
    path "*.log"                             , emit: log
    path "versions.yml"                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    if (params.convert_dotd) {
        extensionconversions = ",.d.gz:.mzML,.d.tar.gz:.mzML,d.tar:.mzML,.d.zip:.mzML,.d:.mzML"
    } else {
        extensionconversions = ",.gz:,.tar.gz:,.tar:,.zip:"
    }
    // .dia files are always passed through without conversion (DIA-NN handles them natively)
    // Compressed .dia files are decompressed but keep the .dia extension
    extensionconversions = "${extensionconversions},.dia.gz:.dia,.dia.tar.gz:.dia,.dia.tar:.dia,.dia.zip:.dia"

    """
    ## -t2 since the one-table format parser is broken in OpenMS2.5
    ## -l for legacy behavior to always add sample columns

    parse_sdrf convert-openms \\
        -t2 -l \\
        --extension_convert raw:mzML$extensionconversions \\
        -s ${sdrf} \\
        $args \\
        2>&1 | tee ${sdrf.baseName}_parsing.log

    mv openms.tsv ${sdrf.baseName}_config.tsv
    mv experimental_design.tsv ${sdrf.baseName}_openms_design.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sdrf-pipelines: \$(parse_sdrf --version 2>&1 | awk -F ' ' '{print \$2}')
    END_VERSIONS
    """
}
