################################# Quality control #################################

rule fastqc:
    input:
        fq1 = expand("{fastqdir}/{{smp}}_{fpostfix}1.fastq.gz", fastqdir=FASTQDIR, fpostfix=FASTQPOSTFIX, smp=SAMPLES),
        fq2 = expand("{fastqdir}/{{smp}}_{fpostfix}2.fastq.gz", fastqdir=FASTQDIR, fpostfix=FASTQPOSTFIX, smp=SAMPLES)
    output:
        fastqc_1 = "1.Quality_control/Sequence_QC/FastQC/{smp}_1_fastqc.html",
        fastqc_2 = "1.Quality_control/sequence_QC/FastQC/{smp}_2_fastqc.html"
    conda:
        "envs/quast.yaml"
    threads:
        SAMPLE_THREADS
    params:
        output = "1.Quality_control/Sequence_QC/FastQC"
    shell:
        """
        fastqc -q -t {threads} -o {params.output} {input.fq1}
        fastqc -q -t {threads} -o {params.output} {input.fq2}
        """


rule summary_report:
    input:
        fastqc_1 = expand("1.Quality_control/Sequence_QC/FastQC/{smp}_1_fastqc.html", smp=SAMPLES),
        fastqc_2 = expand("1.Quality_control/Sequence_QC/FastQC/{smp}_2_fastqc.html", smp=SAMPLES)
    output:
        report = "1.Quality_control/Sequence_QC/MultiQC_report.html"
    params:
        output = "1.Quality_control/Sequence_QC/FastQC"
    shell:
        "multiqc {params.output} --filename {output.report} -p --force"


