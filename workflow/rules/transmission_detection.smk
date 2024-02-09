####### Transmission network inference

rule transmission_clustering:
    input:
        matrix = rules.distance.output.dist_matrix,
        p1 = rules.distance_visualization.output.p1,
        p2 = rules.distance_visualization.output.p2,
        p3 = rules.distance_visualization.output.p3,
        p4 = rules.distance_visualization.output.p4
    output:
        temp(expand("5.Transmission_cluster/{dir}/{outfilename}_cluster_{method}_{cutoff}.csv", dir=METHOD_DIR, outfilename=OUTFILENAME, method=METHOD, cutoff=CUTOFF))
    params:
        output_dir = expand("5.Transmission_cluster/{dir}", dir=METHOD_DIR),
        prefix = OUTFILENAME
    log:
        expand("5.Transmission_cluster/{dir}/genome_clustering.log", dir=METHOD_DIR)
    shell:
        "Rscript --vanilla {SCRIPTS}genome_clustering.R --distance {input.matrix} --date {META_DATA_FILE} "
        "--lambda {LAMBDA} --beta {BETA} --output_dir {params.output_dir} --prefix {params.prefix} "
        "--decimal_date {DECIMAL_DATE} --method {METHOD} --cutoff {CUTOFF} 2> {log}"


rule transmission_network:
    input:
        it = rules.transmission_clustering.output,
        im = rules.distance.output.dist_matrix
    output:
        tc = expand("5.Transmission_cluster/{dir}/transmission_clusters.txt", dir=METHOD_DIR),
        ti = expand("5.Transmission_cluster/{dir}/transmission_detection_run_info.txt", dir=METHOD_DIR)
    params:
        output_dir = expand("5.Transmission_cluster/{dir}", dir=METHOD_DIR)
    log:
        expand("5.Transmission_cluster/{dir}/transmission_detection.log", dir=METHOD_DIR)
    shell:
        "python3 {SCRIPTS}run_transmission_detection.py --cluster {input.it} --distance {input.im} --network {NETWORK} "
        "--output {params.output_dir} --date {META_DATA_FILE} --coord {COORDINATE} --method trans --xmin {CLUSTER_SIZE} 2> {log}"


rule clustering_plot:
    input:
        rules.transmission_network.output.tc
    output:
        p1 = expand("5.Transmission_cluster/{dir}/clustering_visualization.pdf", dir=METHOD_DIR),
        p2 = expand("5.Transmission_cluster/{dir}/clustering_visualization.png", dir=METHOD_DIR),
        p3 = expand("5.Transmission_cluster/{dir}/cluster_member_visualization.pdf", dir=METHOD_DIR),
        p4 = expand("5.Transmission_cluster/{dir}/cluster_member_visualization.png", dir=METHOD_DIR)
    params:
        working_dir = expand("5.Transmission_cluster/{dir}", dir=METHOD_DIR)
    log:
        expand("5.Transmission_cluster/{dir}/clustering_visualization.log", dir=METHOD_DIR)
    shell:
        "Rscript --vanilla {SCRIPTS}cluster_plot.R --input {input} --working_dir {params.working_dir} 2> {log}"

