rule merge_reads:
    input:
        fq1 = "temp/adapt_clip/{smp}/{smp}_1.fastq.gz",
        fq2 = "temp/adapt_clip/{smp}/{smp}_2.fastq.gz"
    output:
        temp("temp/merged/{smp}/{smp}.extendedFrags.fastq.gz"),
        temp("temp/merged/{smp}/{smp}.notCombined_1.fastq.gz"),
        temp("temp/merged/{smp}/{smp}.notCombined_2.fastq.gz"),
        temp("temp/merged/{smp}/{smp}.hist"),
        temp("temp/merged/{smp}/{smp}.histogram")
    params:
        sample = "{smp}",
        overlap = FLASH_OVERLAP,
        outdir = "temp/merged/{smp}",
        threads = SAMPLE_THREADS,
        max_overlap = lambda wildcards: read_length_from_name(wildcards)
    threads:
        SAMPLE_THREADS
    log:
        "temp/merged/{smp}/flash.log"
    shell:
        "flash -t {threads} -z -m {params.overlap} -M {params.max_overlap[0]} -d {params.outdir} "
        "-o {params.sample} {input.fq1} {input.fq2} 2>&1 > {log}"


rule clip:
    input:
        "temp/merged/{smp}/{smp}.notCombined_{nr, \d}.fastq.gz"
    output:
        temp("temp/clipping/{smp}/{smp}.notCombined_clipped_{nr}.fastq")
    log:
        "temp/clipping/{smp}/clip{nr}.log"
    params:
        clip = CLIPPING
    shell:
        "seqtk trimfq -e {params.clip} {input} > {output} 2> {log}"
