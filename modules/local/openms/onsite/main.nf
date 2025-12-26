process ONSITE {
    tag "$meta.mzml_id"
    label 'process_medium'
    label 'openms'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyonsite:0.0.2--pyhdfd78af_0' :
        'quay.io/biocontainers/pyonsite:0.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(mzml_file), path(id_file)

    output:
    tuple val(meta), path("${id_file.baseName}_onsite_*.idXML"), emit: ptm_in_id_onsite
    path "versions.yml", emit: versions
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"

    // Select algorithm (ascore, phosphors, or lucxor)
    def algorithm = params.onsite_algorithm ?: 'lucxor'

    // Basic parameters
    def fragment_tolerance = params.onsite_fragment_tolerance ?: 0.5
    def fragment_units = params.onsite_fragment_error_units ?: 'Da'
    def threads = params.onsite_threads ?: task.cpus
    def add_decoys = params.onsite_add_decoys ?: false
    def min_psms = params.onsite_min_psms ?: 5
    def disable_split_by_charge = params.onsite_disable_split_by_charge ?: false
    def compute_all_scores = params.onsite_compute_all_scores != null ? params.onsite_compute_all_scores : true

    // Algorithm-specific parameters
    def fragment_method = params.onsite_fragment_method ?: meta.dissociationmethod
    def neutral_losses = params.onsite_neutral_losses ? "--neutral-losses ${params.onsite_neutral_losses}" : ""
    def decoy_mass = params.onsite_decoy_mass ? "--decoy-mass ${params.onsite_decoy_mass}" : ""
    def decoy_losses = params.onsite_decoy_neutral_losses ? "--decoy-neutral-losses ${params.onsite_decoy_neutral_losses}" : ""
    def min_psms_param = "--min-num-psms-model ${min_psms}"

    // Debug options - onsite only accepts --debug flag without value
    def debug = params.onsite_debug ? "--debug" : ""

    // Build algorithm-specific parameter strings
    def tolerance_param = ""
    def method_param = ""
    def algorithm_specific_params = ""
    def decoy_param = ""

    if (algorithm == 'lucxor') {
        // LucXor uses --fragment-error-units and --fragment-method
        tolerance_param = "--fragment-error-units ${fragment_units}"
        method_param = fragment_method ? "--fragment-method ${fragment_method}" : ""
        algorithm_specific_params = "${neutral_losses} ${decoy_mass} ${decoy_losses} ${min_psms_param}"
        
        // Add LucXor-specific parameters
        // Note: disable_split_by_charge is only supported by LucXor
        if (disable_split_by_charge) {
            algorithm_specific_params += " --disable-split-by-charge"
        }
        if (compute_all_scores) {
            algorithm_specific_params += " --compute-all-scores"
        }

        // LucXor uses --target-modifications
        // Build target modifications list from params.mod_localization
        if (params.mod_localization) {
            def target_mods = params.mod_localization.tokenize(',').collect { it.trim() }

            // Add decoy modification if enabled
            if (add_decoys) {
                target_mods.add('PhosphoDecoy(A)')
            }

            // Format as command line arguments
            decoy_param = "--target-modifications '${target_mods.join(',')}'"
        } else if (add_decoys) {
            // If no mod_localization specified but decoys enabled, use default with decoy
            decoy_param = "--target-modifications 'Phospho(S),Phospho(T),Phospho(Y),PhosphoDecoy(A)'"
        }
    } else {
        // AScore and PhosphoRS use --fragment-mass-unit
        tolerance_param = "--fragment-mass-unit ${fragment_units}"
        method_param = ""
        algorithm_specific_params = ""
        
        // Add compute_all_scores parameter for AScore and PhosphoRS
        // Note: disable_split_by_charge is LucXor-specific and not used here
        if (compute_all_scores) {
            algorithm_specific_params += " --compute-all-scores"
        }

        // AScore and PhosphoRS use --add-decoys flag
        if (add_decoys) {
            decoy_param = "--add-decoys"
        }
    }

    """
    onsite ${algorithm} \\
        -in ${mzml_file} \\
        -id ${id_file} \\
        -out ${id_file.baseName}_onsite_${algorithm}.idXML \\
        --fragment-mass-tolerance ${fragment_tolerance} \\
        ${tolerance_param} \\
        --threads ${threads} \\
        ${method_param} \\
        ${decoy_param} \\
        ${algorithm_specific_params} \\
        ${debug} \\
        2>&1 | tee ${id_file.baseName}_onsite_${algorithm}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        onsite: \$(onsite --version 2>&1 | grep -oP 'version \\K[0-9.]+' || echo "0.0.1")
    END_VERSIONS
    """
}
