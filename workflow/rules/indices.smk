############### indices
rule fa_index:
    input:
        "{genome}"
    output:
        "{genome}.fai"
    shell:
        "samtools faidx {input}"


rule fasta_dict:
    input:
        "{fastafile}.fasta"
    output:
        "{fastafile}.dict"
    threads:
        SAMPLE_THREADS
    shell:
        "picard -XX:ParallelGCThreads={threads} CreateSequenceDictionary R={input} O={output} QUIET=true VERBOSITY=ERROR"


rule bam_index:
    input:
        "{bamfile}.bam"
    output:
        "{bamfile}.bam.bai"
    threads:
        SAMPLE_THREADS
    shell:
        "samtools index -@ {threads} {input}"


rule bwa_index:
    input:
        "{genome}"
    output:
        "{genome}.bwt",
        "{genome}.amb",
        "{genome}.ann",
        "{genome}.pac",
        "{genome}.sa"
    shell:
        "bwa index {input}"


rule tabix:
    input:
        "{vcffile}.vcf"
    output:
        gz="{vcffile}.vcf.gz",
        tbi="{vcffile}.vcf.gz.tbi"
    shell:
        """
        bgzip -c {input} > {output.gz}
        tabix {output.gz}
        """
