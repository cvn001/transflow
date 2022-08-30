#!/usr/bin/env bash

#####################################################################################################################

display_help() {
    echo "TransFlow - a snakemake pipeline for WGS based MTB transmission detection"
    echo ""
        echo "Usage: bash transflow.sh --configfile configfile.txt [Snakemake options]"
        echo
    echo " Useful Snakemake parameters:"
        echo "    -n,     --dryrun              do not execute anything"
        echo "    -p,     --printshellcmds      print out the shell commands that will be executed"
        echo "    -t,     --timestamp           add a timestamp to all logging output"
        echo "    -c N,   --cores N             use at most N cores in parallel"
        echo "    --ri,   --rerun-incomplete    re-run all jobs the output of which is recognized as incomplete"
        echo "    -q,     --quiet               do not output certain information"
        echo "            --cleanup-shadow      Cleanup old shadow directories which have not been deleted due to failures or power loss"
        echo "            --verbose             print detailed stack traces and detailed operations"
        echo "            --nocolor             Do not use a colored output"
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

mkdir -p "runtime_statistics"

#####################################################################################################################

echo "1st step: Generate quality control statistics (based on fastqc and MultiQC) ..."
\time -v -o "./runtime_statistics/quality_control.txt" snakemake -s "$(dirname "$0")/workflow/quality_control.snakefile" "$@"

#####################################################################################################################

echo "2nd step: variants calling (based on MTB pan-genome & bwa + GATK3) ..."
\time -v -o "./runtime_statistics/variant_calling.txt" snakemake -s "$(dirname "$0")/workflow/variant_calling.snakefile" "$@"

#####################################################################################################################

echo "3rd step: sample clustering and transmission events predicting (based on transcluster + SeqTrack) ..."
\time -v -o "./runtime_statistics/transmission_analysis.txt" snakemake -s "$(dirname "$0")/workflow/transmission_analysis.snakefile" "$@"

#####################################################################################################################

echo "4th step: generating a summary report (html format)..."
\time -v -o "./runtime_statistics/report_generating.txt" snakemake -s "$(dirname "$0")/workflow/report_generating.snakefile" "$@"

#####################################################################################################################
