################### vcf filtering ####################

rule split_vars:
    input:
        vcf = "{vcffile}.vcf",
        dict = GENOME_FDIR + "/" + GENOME + ".dict",
        fai = GENOME_FDIR + "/" + GENOME + ".fasta.fai",
        fasta = GENOME_FDIR + "/" + GENOME + ".fasta"
    output:
        snps = "{vcffile}_snps.vcf",
        indels = "{vcffile}_indels.vcf",
        svs = "{vcffile}_svs.vcf",
        snpsidx = temp("{vcffile}_snps.vcf.idx"),
        indelsidx = temp("{vcffile}_indels.vcf.idx"),
        svsidx = temp("{vcffile}_svs.vcf.idx")
    log:
        "{vcffile}.log"
    shell:
        """
        gatk3 -T SelectVariants -o {output.snps} -V {input.vcf} -selectType SNP -R {input.fasta} &> {log}
        gatk3 -T SelectVariants -o {output.indels} -V {input.vcf} -selectType INDEL -R {input.fasta} &>> {log}
        gatk3 -T SelectVariants -o {output.svs} -V {input.vcf} --selectTypeToExclude INDEL --selectTypeToExclude SNP -R {input.fasta} &>> {log}
        """


## deletions are reported on the previous base with "chr pos . AG A" --> $2 instead of $2-1
rule deletions_bed:
    input:
        "3.SNP_calling/{smp}/{genome}/{smp}_vars_indels.vcf"
    output:
        "3.SNP_calling/{smp}/{genome}/{smp}_vars_deletions.bed"
    shell:
        "grep -v '#' {input} | awk 'BEGIN{{OFS=\"\t\";}} {{if(length($4)>length($5)){{ print $1, $2, $2+(length($4)-length($5)) }}}}' > {output}"


## special filters with python (AF calculation) "chr pos .  A  G" --> $2-1  (vcf is one-based, bed files are zero-based)
rule snps_amb_lowqual_regions:
    input:
        "3.SNP_calling/{smp}/{genome}/{smp}_vars_snps.vcf"
    output:
        "3.SNP_calling/{smp}/{genome}/{smp}_vars_snps_amb-lowqual.bed"
    shell:
        "python3 {SCRIPTS}filter_af_gatk_bed.py -i {input} -o {output} -a {AF_threshold}"


## GATK LowQual filter
rule snps_lowqual_regions:
    input:
        "3.SNP_calling/{smp}/{genome}/{smp}_vars_snps.vcf"
    output:
        "3.SNP_calling/{smp}/{genome}/{smp}_vars_snps_lowqual.bed"
    shell:
        "python3 {SCRIPTS}filter_dp_gatk_bed.py -i {input} -o {output} -d {DP_threshold}"
