################################# Quality control #################################

rule fastqc_before_trimming:
    input:
        fq1 = expand("{fastqdir}/{{smp}}_{fpostfix}1.fastq.gz", fastqdir=FASTQDIR, fpostfix=FASTQPOSTFIX, smp=SAMPLES),
        fq2 = expand("{fastqdir}/{{smp}}_{fpostfix}2.fastq.gz", fastqdir=FASTQDIR, fpostfix=FASTQPOSTFIX, smp=SAMPLES)
    output:
        fastqc_1 = "1.Quality_control/Sequence_QC/Before_trimming/FastQC/{smp}_1_fastqc.html",
        fastqc_2 = "1.Quality_control/Sequence_QC/Before_trimming/FastQC/{smp}_2_fastqc.html"
    threads:
        SAMPLE_THREADS
    params:
        output = "1.Quality_control/Sequence_QC/Before_trimming/FastQC"
    shell:
        """
        fastqc -q -t {threads} -o {params.output} {input.fq1}
        fastqc -q -t {threads} -o {params.output} {input.fq2}
        """


rule summary_report_before:
    input:
        fastqc_1 = expand("1.Quality_control/Sequence_QC/Before_trimming/FastQC/{smp}_1_fastqc.html", smp=SAMPLES),
        fastqc_2 = expand("1.Quality_control/Sequence_QC/Before_trimming/FastQC/{smp}_2_fastqc.html", smp=SAMPLES)
    output:
        report = "1.Quality_control/Sequence_QC/Before_trimming/MultiQC_report.html"
    params:
        output = "1.Quality_control/Sequence_QC/Before_trimming/FastQC"
    shell:
        "multiqc {params.output} -n {output.report} -p -q --force"


rule fastqc_after_trimming:
    input:
        fq1 = "temp/ready/{smp}/{smp}.pe_1.fastq.gz",
        fq2 = "temp/ready/{smp}/{smp}.pe_2.fastq.gz",
        fqm = "temp/ready/{smp}/{smp}.se.fastq.gz"
    output:
        fastqc_1 = "1.Quality_control/Sequence_QC/After_trimming/FastQC/{smp}.pe_1_fastqc.html",
        fastqc_2 = "1.Quality_control/Sequence_QC/After_trimming/FastQC/{smp}.pe_2_fastqc.html",
        fastqc_m = "1.Quality_control/Sequence_QC/After_trimming/FastQC/{smp}.se_fastqc.html"
    threads:
        SAMPLE_THREADS
    params:
        output = "1.Quality_control/Sequence_QC/After_trimming/FastQC"
    shell:
        """
        fastqc -q -t {threads} -o {params.output} {input.fq1}
        fastqc -q -t {threads} -o {params.output} {input.fq2}
        fastqc -q -t {threads} -o {params.output} {input.fqm}
        """


rule summary_report_after:
    input:
        fastqc_1 = expand("1.Quality_control/Sequence_QC/After_trimming/FastQC/{smp}.pe_1_fastqc.html", smp=SAMPLES),
        fastqc_2 = expand("1.Quality_control/Sequence_QC/After_trimming/FastQC/{smp}.pe_2_fastqc.html", smp=SAMPLES),
        fastqc_m = expand("1.Quality_control/Sequence_QC/After_trimming/FastQC/{smp}.se_fastqc.html", smp=SAMPLES)
    output:
        report = "1.Quality_control/Sequence_QC/After_trimming/MultiQC_report.html"
    params:
        output = "1.Quality_control/Sequence_QC/After_trimming/FastQC"
    shell:
        "multiqc {params.output} -n {output.report} -p -q --force"
