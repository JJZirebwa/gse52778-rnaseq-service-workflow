nextflow.enable.dsl = 2


process VALIDATE_SERVICE_INTAKE {
    tag 'GSE52778 intake'

    input:
    path samplesheet
    path intake_validator

    output:
    tuple path('validated_samples.csv'), path('intake_validation.json')

    script:
    """
    cp ${samplesheet} validated_samples.csv
    python3 ${intake_validator} \\
      --samplesheet validated_samples.csv \\
      --output intake_validation.json
    """
}


process RUN_DESEQ2_PAIRED {
    tag 'paired ~ cell + dex'
    container params.container_image

    input:
    tuple path(samplesheet), path(intake_validation)
    path analysis_script

    output:
    tuple \
        path('validated_samples.csv'), \
        path('intake_validation.json'), \
        path('full_results.csv'), \
        path('fdr_0.05_results.csv'), \
        path('analysis_metadata.json'), \
        path('library_sizes.png'), \
        path('count_distributions.png'), \
        path('pca.png'), \
        path('sample_distance_heatmap.png'), \
        path('sessionInfo.txt')

    script:
    """
    Rscript --vanilla ${analysis_script} \\
      --samplesheet validated_samples.csv \\
      --output-dir .
    """
}


process VERIFY_RESULTS {
    tag 'DE result contract'

    input:
    tuple \
        path(samplesheet), \
        path(intake_validation), \
        path(full_results), \
        path(fdr_results), \
        path(analysis_metadata), \
        path(library_sizes), \
        path(count_distributions), \
        path(pca), \
        path(sample_distance), \
        path(session_info)
    path result_validator

    output:
    tuple \
        path('validated_samples.csv'), \
        path('intake_validation.json'), \
        path('full_results.csv'), \
        path('fdr_0.05_results.csv'), \
        path('analysis_metadata.json'), \
        path('library_sizes.png'), \
        path('count_distributions.png'), \
        path('pca.png'), \
        path('sample_distance_heatmap.png'), \
        path('sessionInfo.txt'), \
        path('result_validation.json')

    script:
    """
    python3 ${result_validator} \\
      --full-table full_results.csv \\
      --fdr-table fdr_0.05_results.csv \\
      --metadata analysis_metadata.json \\
      --output result_validation.json
    """
}


process BUILD_HANDOFF_BUNDLE {
    tag 'service handoff'
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    tuple \
        path(samplesheet), \
        path(intake_validation), \
        path(full_results), \
        path(fdr_results), \
        path(analysis_metadata), \
        path(library_sizes), \
        path(count_distributions), \
        path(pca), \
        path(sample_distance), \
        path(session_info), \
        path(result_validation)
    path handover_document
    path interpretation_document
    path bundle_builder

    output:
    path 'handoff'

    script:
    """
    python3 ${bundle_builder} \\
      --output-dir handoff \\
      --git-commit '${params.git_commit}' \\
      --container-image '${params.container_image}' \\
      --base-container '${params.base_container}' \\
      --files \\
        ${samplesheet} \\
        ${intake_validation} \\
        ${full_results} \\
        ${fdr_results} \\
        ${analysis_metadata} \\
        ${library_sizes} \\
        ${count_distributions} \\
        ${pca} \\
        ${sample_distance} \\
        ${session_info} \\
        ${result_validation} \\
        ${handover_document} \\
        ${interpretation_document}
    """
}


workflow {
    samplesheet = file(params.samplesheet, checkIfExists: true)
    intake_validator = file("${projectDir}/scripts/validate_service_intake.py", checkIfExists: true)
    result_validator = file("${projectDir}/scripts/verify_results.py", checkIfExists: true)
    analysis_script = file("${projectDir}/scripts/run_deseq2.R", checkIfExists: true)
    bundle_builder = file("${projectDir}/scripts/build_handoff.py", checkIfExists: true)
    handover_document = file("${projectDir}/docs/HANDOVER.md", checkIfExists: true)
    interpretation_document = file("${projectDir}/docs/INTERPRETATION.md", checkIfExists: true)

    intake = VALIDATE_SERVICE_INTAKE(samplesheet, intake_validator)
    analysis = RUN_DESEQ2_PAIRED(intake, analysis_script)
    verified = VERIFY_RESULTS(analysis, result_validator)
    BUILD_HANDOFF_BUNDLE(verified, handover_document, interpretation_document, bundle_builder)
}
