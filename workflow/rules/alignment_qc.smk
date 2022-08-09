rule qualimap:
    input:
        expand("temp/ready_forvar/{{smp}}/{genome}/{{smp}}_ready.bam", smp=SAMPLES, genome=GENOME)
    output:
        d = directory("1.Quality_control/Alignment_QC/Qualimap/{smp}"),
        f = "1.Quality_control/Alignment_QC/Qualimap/{smp}/qualimapReport.html"
    # threads:
    #     SAMPLE_THREADS
    log:
        "1.Quality_control/Alignment_QC/Qualimap/{smp}/qualimap.log"
    shell:
        "qualimap bamqc -bam {input} -nt 1 -outdir {output.d} &> {log}"


rule alignment_summary:
    input:
        q = expand("1.Quality_control/Alignment_QC/Qualimap/{smp}/qualimapReport.html", smp=SAMPLES)
    output:
        report = "1.Quality_control/Alignment_QC/MultiQC_report.html"
    params:
        output = "1.Quality_control/Alignment_QC/Qualimap",
        cl = "qualimap_config: { general_stats_coverage: [1,10] }"
    shell:
        'multiqc {params.output} -n {output.report} --force -p -q --cl_config "{params.cl}"'
