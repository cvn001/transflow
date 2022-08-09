##### MTBC identification and reads filtering

rule kraken_detect:
    input:
        fq1 = rules.adapter_clip.output.fq1,
        fq2 = rules.adapter_clip.output.fq2
    output:
        reads = temp("2.MTBC_identification/{smp}/{smp}.kraken"),
        report = "2.MTBC_identification/{smp}/{smp}.report"
    log:
        "2.MTBC_identification/{smp}/{smp}.log"
    threads:
        SAMPLE_THREADS
    params:
        kraken_db = config["kraken_db"]
    shell:
        """
        kraken --db {params.kraken_db} --threads {threads} --gzip-compressed --fastq-input --paired {input.fq1} {input.fq2} > {output.reads} 2> {log}
        kraken-report --db {params.kraken_db} {output.reads} > {output.report} 2>> {log}
        """

if MTBC_READS_ONLY:
    rule taxonomy_translate:
        input:
            reads = rules.kraken_detect.output.reads
        output:
            labels = temp("2.MTBC_identification/{smp}/{smp}.labels")
        log:
            "2.MTBC_identification/{smp}/{smp}.log"
        params:
            kraken_db = config["kraken_db"]
        shell:
            "kraken-translate {input.reads} --db {params.kraken_db} > {output.labels} 2>> {log}"


    rule make_read_list:
        input:
            rules.taxonomy_translate.output.labels
        output:
            read_list = temp("2.MTBC_identification/{smp}/{smp}.filtered.readlist")
        log:
            "2.MTBC_identification/{smp}/{smp}.log"
        params:
            taxon = "Mycobacterium tuberculosis complex"
        shell:
            "fgrep '{params.taxon}' {input} | cut -f1 > {output} 2>> {log}"


    rule pick_reads:
        input:
            fq1 = rules.adapter_clip.output.fq1,
            fq2 = rules.adapter_clip.output.fq2,
            read_list = rules.make_read_list.output.read_list
        output:
            flt_fq1 = temp("temp/adapt_clip/{smp}/{smp}_cleaning_mtbc_1.fastq"),
            flt_fq2 = temp("temp/adapt_clip/{smp}/{smp}_cleaning_mtbc_2.fastq")
        log:
            "2.MTBC_identification/{smp}/{smp}.log"
        shell:
            """
            seqtk subseq {input.fq1} {input.read_list} > {output.flt_fq1} 2>> {log}
            seqtk subseq {input.fq2} {input.read_list} > {output.flt_fq2} 2>> {log}
            """
    
    rule preparing_reads:
        input:
            fq1 = rules.pick_reads.output.flt_fq1,
            fq2 = rules.pick_reads.output.flt_fq2
        output:
            fq1_gz = temp("temp/adapt_clip/{smp}/{smp}_cleaning_mtbc_1.fastq.gz"),
            fq2_gz = temp("temp/adapt_clip/{smp}/{smp}_cleaning_mtbc_2.fastq.gz")
        log:
            "temp/adapt_clip/{smp}/preparing_mtbc_reads.log"
        shell:
            """
            gzip -c {input.fq1} > {output.fq1_gz}
            gzip -c {input.fq2} > {output.fq2_gz}
            """


rule kraken_filter:
    input:
        rules.kraken_detect.output.report
    output:
        "2.MTBC_identification/{smp}/{smp}.result"
    params:
        sample_name = "{smp}"
    shell:
        "python {SCRIPTS}kraken_report_extract.py {input} {output} {params.sample_name}"


rule kraken_filter_merge:
    input:
        expand("2.MTBC_identification/{smp}/{smp}.result", smp=SAMPLES)
    output:
        "2.MTBC_identification/kraken.result"
    log:
        "2.MTBC_identification/kraken_result_merge.log"
    shell:
        """
        echo -e "##The percentage of reads mapping to MTBC and the category detected by kraken" > {output}
        echo -e "##Kraken QC threshold: {KRAKEN_FILTER}% to MTBC (samples with lower percentage won't run the rest of workflow)" >> {output}
        echo -e "#sample\tpercentage_MTBC\tcategory_detected" >> {output}
        cat {input} >> {output}
        """


rule mtbc_filter:
    input:
        rules.kraken_filter_merge.output
    output:
        "2.MTBC_identification/MTBC_samples.txt"
    run:
        in_file = "2.MTBC_identification/kraken.result"
        out_file = "2.MTBC_identification/MTBC_samples.txt"
        with open(in_file, 'r') as f, open(out_file, 'w') as o:
            for each_line in f.readlines():
                if '#' not in each_line:
                    s_list = each_line.strip().split('\t')
                    sample = s_list[0]
                    percent = float(s_list[1])
                    if percent >= KRAKEN_FILTER:
                        o.write(sample + '\n')
                    else:
                        print("Warning! [{}] is not a MTBC strain.".format(sample))
