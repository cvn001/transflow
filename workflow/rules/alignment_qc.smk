rule qualimap:
    input:
        expand("temp/ready_forvar/{{smp}}/{genome}/{{smp}}_ready.bam", smp=SAMPLES, genome=GENOME)
    output:
        d = directory("1.Quality_control/Alignment_QC/Qualimap/{smp}"),
        f = "1.Quality_control/Alignment_QC/Qualimap/{smp}/qualimapReport.html"
    threads:
        SAMPLE_THREADS
    shell:
        "qualimap bamqc -bam {input} -nt {threads} --java-mem-size=2G -outdir {output.d}"


rule alignment_summary:
    input:
        q = expand("1.Quality_control/Alignment_QC/Qualimap/{smp}/qualimapReport.html", smp=SAMPLES)
    output:
        report = "1.Quality_control/Alignment_QC/MultiQC_report.html"
    params:
        output = "1.Quality_control/Alignment_QC/Qualimap",
        cl = "qualimap_config: { general_stats_coverage: [1,10] }"
    shell:
        'multiqc {params.output} --filename {output.report} --force -p --cl_config "{params.cl}"'
