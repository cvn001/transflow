##### Summary_report #####

rule generate_yaml:
    input:
        WORK_DIR
    output:
        temp(expand('7.Summary_report/{dir}/report.yaml', dir=METHOD_DIR))
    params:
        b = KRAKEN_CUTOFF,
        c = MQ_threshold,
        d = AF_threshold,
        e = DP_threshold,
        f = METHOD,
        g = SNP_CUTOFF,
        i = TRANS_CUTOFF,
        j = METHOD_DIR
    log:
        expand('7.Summary_report/{dir}/generate_yaml.log', dir=METHOD_DIR)
    shell:
        "python3 {SCRIPTS}generate_report_yaml.py -a {input} -b {params.b} -c {params.c} -d {params.d}"
        " -e {params.e} -f {params.f} -g {params.g} -i {params.i} -j {params.j} 2> {log}"

# Generate summary report using a R markdown script
rule generate_report:
    input:
        rules.generate_yaml.output
    output:
        expand("7.Summary_report/{dir}/summary_report.html", dir=METHOD_DIR)
    log:
        expand("7.Summary_report/{dir}/summary_report.log", dir=METHOD_DIR)
    shell:
        # Run a R script to pass parameters to Rmd script
        "Rscript --vanilla {SCRIPTS}run_rmd.R {SCRIPTS}generate_summary_html.Rmd {output} {input} 2> {log}"
        