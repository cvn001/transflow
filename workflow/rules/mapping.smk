############## mapping

rule map_pe:
    input:
        fq1 = "temp/ready/{smp}/{smp}.pe_1.fastq.gz",
        fq2 = "temp/ready/{smp}/{smp}.pe_2.fastq.gz",
        index = GENOME_FDIR + "/{genome}.fasta.bwt",
        g = GENOME_FDIR + "/{genome}.fasta"
    output:
        f = "temp/mapped/{smp}/{genome}/{smp}_pe_Aligned.out.bam"
    log:
        "temp/mapped/{smp}/{genome}/{smp}_pe_bwa_mem.log"
    threads:
        SAMPLE_THREADS
    shell:
        "bwa mem -v 1 -t {threads} {input.g} {input.fq1} {input.fq2} 2> {log} | samtools view -Sb - > {output}"


rule map_se:
    input:
        fq = "temp/ready/{smp}/{smp}.se.fastq.gz",
        index = GENOME_FDIR + "/{genome}.fasta.bwt",
        g = GENOME_FDIR + "/{genome}.fasta"
    output:
        f = "temp/mapped/{smp}/{genome}/{smp}_se_Aligned.out.bam"
    log:
        "temp/mapped/{smp}/{genome}/{smp}_se_bwa_mem.log"
    threads:
        SAMPLE_THREADS
    shell:
        "bwa mem -v 1 -t {threads} {input.g} {input.fq} 2> {log} | samtools view -Sb - > {output}"
