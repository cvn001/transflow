############### bam file modifications

rule sort_bam:
    input:
        "temp/mapped/{smp}/{genome}/{smp}_{reads}_Aligned.out.bam"
    output:
        "temp/mapped/{smp}/{genome}/{smp}_{reads}_Aligned.sortedByCoord.out.bam"
    threads:
        SAMPLE_THREADS
    log:
        "temp/mapped/{smp}/{genome}/{smp}_{reads}_Aligned.sortedByCoord.out.tmp"
    shell:
        "samtools sort -@ {threads} -T {log} {input} > {output}"


rule combine_bam:
    input:
        se = "temp/mapped/{smp}/{genome}/{smp}_se_Aligned.sortedByCoord.out.bam",
        pe = "temp/mapped/{smp}/{genome}/{smp}_pe_Aligned.sortedByCoord.out.bam"
    output:
        "temp/ready_forvar/{smp}/{genome}/{smp}.bam"
    threads:
        SAMPLE_THREADS
    shell:
        "samtools merge -@ {threads} {output} {input.se} {input.pe}"


rule mark_duplicates:
    input:
        "temp/ready_forvar/{smp}/{genome}/{smp}.bam"
    output:
        bam = "temp/ready_forvar/{smp}/{genome}/{smp}_markeddup.bam",
        txt = "temp/ready_forvar/{smp}/{genome}/{smp}_metrics.txt"
    log:
        "temp/ready_forvar/{smp}/{genome}/mark_duplicates.log"
    threads:
        SAMPLE_THREADS
    shell:
        "picard -XX:ParallelGCThreads={threads} MarkDuplicates I={input} O={output.bam} METRICS_FILE={output.txt} 2>&1 >{log}"


rule add_readgroups:
    input:
        "temp/ready_forvar/{smp}/{genome}/{smp}_markeddup.bam"
    output:
        "temp/ready_forvar/{smp}/{genome}/{smp}_ready.bam"
    log:
        "temp/ready_forvar/{smp}/{genome}/add_readgroups.log"
    threads:
        SAMPLE_THREADS
    run:
        sm = wildcards["smp"]
        lb = "lb"
        pu = "pu"
        shell("picard -XX:ParallelGCThreads={threads} AddOrReplaceReadGroups I={input[0]} O={output} LB={lb} PL=illumina PU={pu} SM={sm} 2>&1 >{log}")
