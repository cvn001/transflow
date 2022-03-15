#!/usr/bin/env bash

display_help() {
    echo "transmission-pipeline - NGS Alignment and SNP Calling pipeline for MTB transmission detection (powered by Snakemake)"
    echo ""
        echo "Usage: transmission-pipeline --configfile configfile.txt [Snakemake options]"
        echo
    echo " Useful Snakemake parameters:"
        echo "   -n,   --dryrun           do not execute anything"
        echo "   -p,   --printshellcmds   print out the shell commands that will be executed"
        echo "   -t,   --timestamp        add a timestamp to all logging output"
        echo "   -j N, --jobs N           use at most N cores in parallel"
    echo
        echo " Full list of parameters:"
    echo "   --help                 show Snakemake help (or snakemake -h)"
        echo
    exit 0
}

# shellcheck disable=SC2166
if [ "$1" == "" -o "$1" == "-h" -o \( "$1" != "--configfile" -a "$1" != "--help" \) ]; then
  display_help
  exit 0
fi
#####################################################################################################################
echo "1st step: Generate quality control statistics (based on fastqc and MultiQC) ..."
snakemake -s "${BASH_SOURCE%/*}/workflow/quality_control.snakefile" "$@"
#####################################################################################################################
echo "2nd step: MTBC filtering (based on Kraken) ..."
snakemake -s "${BASH_SOURCE%/*}/workflow/mtbc_filtering.snakefile" "$@"
#####################################################################################################################
echo "3rd step: variants calling (based on MTB pan-genome & bwa + GATK3) ..."
snakemake -s "${BASH_SOURCE%/*}/workflow/variant_calling.snakefile" "$@"
#####################################################################################################################
echo "4th step: sample clustering and transmission events predicting (based on transcluster + SeqTrack) ..."
snakemake -s "${BASH_SOURCE%/*}/workflow/transmission_detection.snakefile" "$@"
#####################################################################################################################
echo "5th step: generating a summary report (html format)..."
snakemake -s "${BASH_SOURCE%/*}/workflow/generating_report.snakefile" "$@"
