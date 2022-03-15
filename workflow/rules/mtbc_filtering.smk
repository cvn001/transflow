##### MTBC filtering

rule kraken_detect:
    input:
        fq1 = expand("{fastqdir}/{{smp}}_{fpostfix}1.fastq.gz", fastqdir=FASTQDIR, fpostfix=FASTQPOSTFIX, smp=SAMPLES),
        fq2 = expand("{fastqdir}/{{smp}}_{fpostfix}2.fastq.gz", fastqdir=FASTQDIR, fpostfix=FASTQPOSTFIX, smp=SAMPLES)
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
        kraken --db {params.kraken_db} --threads {threads} --gzip-compressed --fastq-input \
        --paired {input.fq1} {input.fq2} > {output.reads} 2> {log}
        kraken-report --db {params.kraken_db} {output.reads} > {output.report} 2>> {log}  
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
