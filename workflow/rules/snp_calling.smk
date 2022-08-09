############### snp calling

rule bam_filter:
    input:
        "{bamfile}.bam"
    output:
        temp("{bamfile}_mqfiltered.bam")
    shell:
        "samtools view -q {MQ_threshold} -b {input} > {output}"

rule snp_calling:
    input:
        fai = GENOME_FDIR + "/{genome}.fasta.fai",
        fasta = GENOME_FDIR + "/{genome}.fasta",
        dict = GENOME_FDIR + "/{genome}.dict",
        bam = "temp/ready_forvar/{smp}/{genome}/{smp}_ready.bam",
        bai = "temp/ready_forvar/{smp}/{genome}/{smp}_ready.bam.bai"
    output:
        vcf = "3.SNP_calling/{smp}/{genome}/{smp}.g.vcf",
        idx = "3.SNP_calling/{smp}/{genome}/{smp}.g.vcf.idx"
    threads:
        SAMPLE_THREADS
    log:
        "3.SNP_calling/{smp}/{genome}/snp_calling.log"
    shell:
        "gatk3 -T HaplotypeCaller -o {output.vcf} -I {input.bam} --sample_ploidy 2 -nct {threads} -R {input.fasta} -mmq {MQ_threshold} -ERC GVCF &> {log}"


rule genotype_calling:
    input:
        gvcf = "3.SNP_calling/{smp}/{genome}/{smp}.g.vcf",
        fai = GENOME_FDIR + "/{genome}.fasta.fai",
        fasta = GENOME_FDIR + "/{genome}.fasta",
        dict = GENOME_FDIR + "/{genome}.dict"
    output:
        vcf = "3.SNP_calling/{smp}/{genome}/{smp}.vcf",
        idx = "3.SNP_calling/{smp}/{genome}/{smp}.vcf.idx"
    threads:
        SAMPLE_THREADS
    log:
        "3.SNP_calling/{smp}/{genome}/genotyping.log"
    shell:
        "gatk3 -T GenotypeGVCFs -R {input.fasta} --variant {input.gvcf} -o {output.vcf} -nt {threads} --sample_ploidy 2 --includeNonVariantSites -stand_call_conf {DP_threshold} &> {log}"


rule gatk_vars:
    input:
        "3.SNP_calling/{smp}/{genome}/{smp}.vcf"
    output:
        vcf = "3.SNP_calling/{smp}/{genome}/{smp}_vars.vcf",
        bed = "3.SNP_calling/{smp}/{genome}/{smp}_uncalled_singlepos.bed",
    shell:
        "python3 {SCRIPTS}filter_gatk.py -i {input} --output_vcf {output.vcf} --output_bed {output.bed}"


rule uncalled:
    input:
        "3.SNP_calling/{smp}/{genome}/{smp}_uncalled_singlepos.bed"
    output:
        "3.SNP_calling/{smp}/{genome}/{smp}_uncalled.bed"
    shell:
        "bedtools merge -i {input} > {output}"


rule genome_coverage:
    input:
        bam = "temp/ready_forvar/{smp}/{genome}/{smp}_ready_mqfiltered.bam",
        fasta = GENOME_FDIR + "/{genome}.fasta"
    output:
        temp("3.SNP_calling/{smp}/{genome}/{smp}_coverage.bed")
    shell:
        "bedtools genomecov -bga -ibam {input.bam} > {output}"


rule low_coverage:
    input:
        "3.SNP_calling/{smp}/{genome}/{smp}_coverage.bed"
    output:
        "3.SNP_calling/{smp}/{genome}/{smp}_uncovered.bed"
    shell:
        "awk '$4 < {DP_threshold}' {input} | bedtools merge -i - > {output}"
