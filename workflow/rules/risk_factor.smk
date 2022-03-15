############### Transmission risk factor testing


rule risk_factor_testing:
    input:
        a = expand("5.Transmission_cluster/{output}/transmission_clusters.txt", output=METHOD_DIR),
        b = expand("5.Transmission_cluster/{output}/transmission_detection_run_info.txt", output=METHOD_DIR),
        c = expand("5.Transmission_cluster/{output}/clustering_visualization.pdf", output=METHOD_DIR),
        d = expand("5.Transmission_cluster/{output}/clustering_visualization.png", output=METHOD_DIR),
        e = expand("5.Transmission_cluster/{output}/cluster_member_visualization.pdf", output=METHOD_DIR),
        f = expand("5.Transmission_cluster/{output}/cluster_member_visualization.png", output=METHOD_DIR)
    output:
        expand("6.Risk_factor/{output}/univariable_regression_result.html", output=METHOD_DIR),
        expand("6.Risk_factor/{output}/univariable_regression_result.rtf", output=METHOD_DIR),
        expand("6.Risk_factor/{output}/characteristic_data.txt", output=METHOD_DIR)
    params:
        workdir = expand("6.Risk_factor/{output}", output=METHOD_DIR)
    log:
        expand("6.Risk_factor/{output}/risk_factor_testing.log", output=METHOD_DIR)
    shell:
        "Rscript --vanilla {SCRIPTS}transmission_risk_factor.R -k {META_DATA_FILE} -c {input.a}"
        " -l {CHARACTERS} -o {params.workdir} 2> {log}"
