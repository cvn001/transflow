rule getinfo:
    input:
         snp = "3.SNP_calling",
         gaps = GENOMEGAPSFILE,
         ids = SAMPLE_LIST,
         fai = expand("{genomefile}.fai", genomefile=GENOME_FILE)
    output:
         expand("4.SNP_distance/{outfilename}.RSave", outfilename=OUTFILENAME)
    shell:
        "Rscript --vanilla {SCRIPTS}get_snp_summary.R --input_dir {input.snp} --id_file {input.ids} "
        "--fastaidx_file {input.fai} --genome_gaps_file {input.gaps} --output_file {output}"


rule distance:
    input:
        dist = rules.getinfo.output
    output:
        dist_matrix = expand("4.SNP_distance/{outfilename}_pairwise_distance_matrix.txt", outfilename=OUTFILENAME),
        snp_table = expand("4.SNP_distance/{outfilename}_snp_matrix.txt", outfilename=OUTFILENAME),
        pair_count = temp(expand("4.SNP_distance/samples_pairwise_count.txt", outfilename=OUTFILENAME)),
    params:
        prefix = "samples",
        exclude = check_exclude,
        output = "4.SNP_distance"
    log:
        "4.SNP_distance/distance.log"
    shell:
        "Rscript --vanilla {SCRIPTS}analyze_snp_summary_matrices_pairwise.R --input_data {input.dist} "
        "--output_dir {params.output} --output_prefix {params.prefix} {params.exclude} 2> {log}"


rule distance_visualization:
    input:
        rules.distance.output.dist_matrix
    output:
        p1 = expand("4.SNP_distance/{outfilename}_pairwise_distance_distribution.pdf", outfilename=OUTFILENAME),
        p2 = expand("4.SNP_distance/{outfilename}_pairwise_distance_distribution.png", outfilename=OUTFILENAME),
        p3 = expand("4.SNP_distance/{outfilename}_pairwise_distance_heatmap.pdf", outfilename=OUTFILENAME),
        p4 = expand("4.SNP_distance/{outfilename}_pairwise_distance_heatmap.png", outfilename=OUTFILENAME)
    params:
        prefix = "samples",
        output = "4.SNP_distance",
        show_rownames = "True",
        show_colnames = "True"
    log:
        "4.SNP_distance/visualization.log"
    shell:
        "Rscript --vanilla {SCRIPTS}distance_visualization.R {input} {params.output} {params.show_colnames} "
        "{params.show_rownames} {params.prefix} 2> {log}"
