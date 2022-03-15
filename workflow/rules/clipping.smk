############## clipping

rule adapter_clip:
    input:
        fq1 = expand("{fastqdir}/{{smp}}_{fpostfix}1.fastq.gz", fastqdir=FASTQDIR, fpostfix=FASTQPOSTFIX, smp=SAMPLES),
        fq2 = expand("{fastqdir}/{{smp}}_{fpostfix}2.fastq.gz", fastqdir=FASTQDIR, fpostfix=FASTQPOSTFIX, smp=SAMPLES),
        af = ADAPTER_FILE
    output:
        fq1 = temp("temp/adapt_clip/{smp}/{smp}_1.fastq.gz"),
        fq2 = temp("temp/adapt_clip/{smp}/{smp}_2.fastq.gz"),
        fqu1 = temp("temp/adapt_clip/{smp}/{smp}_1_up.fastq.gz"),
        fqu2 = temp("temp/adapt_clip/{smp}/{smp}_2_up.fastq.gz")
    log:
        "temp/adapt_clip/{smp}/illuminaclip.log"
    threads:
        SAMPLE_THREADS
    shell:
        "trimmomatic PE -phred33 -quiet -threads {threads} {input.fq1} {input.fq2} {output.fq1} {output.fqu1} {output.fq2} {output.fqu2} ILLUMINACLIP:{input.af}:2:30:10:1:true 2> {log}"
