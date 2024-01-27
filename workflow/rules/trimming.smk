############# NGS data trimming

rule qc_pe:
    input:
        notcomb1 = "temp/clipping/{smp}/{smp}.notCombined_clipped_1.fastq",
        notcomb2 = "temp/clipping/{smp}/{smp}.notCombined_clipped_2.fastq"
    output:
        p1 = temp("temp/QCPE/{smp}/r1_P_{smp}.notCombined_clipped_1.fastq"),
        up1 = temp("temp/QCPE/{smp}/r1_UP_{smp}.notCombined_clipped_1.fastq"),
        p2 = temp("temp/QCPE/{smp}/r2_P_{smp}.notCombined_clipped_2.fastq"),
        up2 = temp("temp/QCPE/{smp}/r2_UP_{smp}.notCombined_clipped_2.fastq")
    log:
        "temp/QCPE/{smp}/QCPE.log"
    threads:
        SAMPLE_THREADS
    shell:
        "trimmomatic PE -quiet -threads {threads} -phred33 {input.notcomb1} {input.notcomb2} {output.p1} {output.up1} {output.p2} {output.up2} SLIDINGWINDOW:4:{MIN_QUAL} MINLEN:{READ_MIN_LEN} 2>&1 > {log}"


rule qc_se:
    input:
        extfraq = "temp/merged/{smp}/{smp}.extendedFrags.fastq.gz"
    output:
        temp("temp/QCSE/{smp}/{smp}.extendedFrags.fastq")
    log:
        "temp/QCSE/{smp}/QCSE.log"
    threads:
        SAMPLE_THREADS
    shell:
        "trimmomatic SE -quiet -threads {threads} -phred33 {input.extfraq} {output} SLIDINGWINDOW:4:{MIN_QUAL} MINLEN:{READ_MIN_LEN} 2>&1 > {log}"


rule prepare_mapping:
    input:
        fq1up = "temp/QCPE/{smp}/r1_UP_{smp}.notCombined_clipped_1.fastq",
        fq2up = "temp/QCPE/{smp}/r2_UP_{smp}.notCombined_clipped_2.fastq",
        fq1p = "temp/QCPE/{smp}/r1_P_{smp}.notCombined_clipped_1.fastq",
        fq2p = "temp/QCPE/{smp}/r2_P_{smp}.notCombined_clipped_2.fastq",
        fqm = "temp/QCSE/{smp}/{smp}.extendedFrags.fastq"
    output:
        fq1 = "temp/ready/{smp}/{smp}.pe_1.fastq.gz",
        fq2 = "temp/ready/{smp}/{smp}.pe_2.fastq.gz",
        fqm = "temp/ready/{smp}/{smp}.se.fastq.gz"
    shell:
        """
        gzip -c {input.fq1p} > {output.fq1}
        gzip -c {input.fq2p} > {output.fq2}
        cat {input.fqm} {input.fq1up} {input.fq2up} | gzip -c > {output.fqm}
        """
