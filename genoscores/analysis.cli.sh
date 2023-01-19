#!/usr/bin/env bash

#' This script computes scores for specified intermediate trait batching target
#' genotypes by chromosome.
#'
#' Intermediate traits useful to test omnigenic hypothesis:
#' eQTLs (eQTLGen); studyid=36,
#' pQTLs (Somalogic 2018); studyid=26,
#' pancreatic islets eQTLs; studyid=1007,
#' immune cell female twins; studyid=35.
## ============================================================================

## Set the path to the output directory to store the GENOSCORES results.
## The outputs for each chromosome will be saved in a respective subdirectory.
## Examples: outputdir="." or outputdir="/home/user/scores"
outputdir=$1

echo "GENOSCORES results will be saved to: $outputdir."

## Set the base name (with full path) of the input file (excluding ID).
## Example: if files are called /home/user/cohort_chr<1-22>.bim|bed|fam,
## set basename="/home/user/cohort_chr"
basename=$2

echo "Target genotypes are in: $basename."

## Set path to the GENOSCORES analysis R script. Here you can use
## example.analysis.R
## Example: gscript="/home/user/example.analysis.R"
gscript=$3

echo "GENOSCORES analysis script: $gscript will be executed."

## Specify the number of processes to be run in parallel, each chromosome will
## be analysed in one process
N=2

## Function to process one chromosome at a time.
## If you target file has suffixes after chromosome ID, for example:
## 'cohort_chr21_filtered', they need to be manually added to inputfile:
## inputfile="${basename}/cohort_chr${chr}_filtered"
processchrom() {
    inputfile="${basename}/cohort_chr${chr}"
    chromdir="${outputdir}/${chr}"

    if [[ ! -d "${chromdir}" ]]; then
        mkdir -p "${chromdir}"
    fi

    logfile="${chromdir}"/log.Rout

    echo "Processing chromosome: ${chr}"
    echo "Input file: ${inputfile}"
    echo "Saving results to: ${chromdir} ..."

    ## Run main GENOSCORES script with CLI options:
    ## chromdir: subdirectory in the output directory to store results from
    ##           the analysis for the current chromosome.
    ## inputfile: PLINK file with the current chromosome.
    ## logfile: file to store the logs from GENOSCORES.
    Rscript --no-save --no-restore --verbose "$(realpath "${gscript}")" \
         "${chromdir}" "${inputfile}" 2> "${logfile}" > /dev/null
    echo "Done."
}

## Run analysis changing chromosome in a given range
for chr in $(seq 1 22); do
    processchrom &

    ## allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        ## now there are $N jobs already running, so wait here for any job
        ## to be finished so there is a place to start next one.
        wait -n
    fi
done

## No more jobs to be started but wait for pending jobs
wait

echo "All analyses have been executed."
