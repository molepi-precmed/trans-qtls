#!/usr/bin/env bash

#' This script is used to run GENOSCORES on a single chromosome in an HPC
#' environment with qsub
## ============================================================================

## Read the path for the directory where the results are to be stored. Typically
## it is on scratch space on HPC
outputdir=$1

## Read the base name (with full path) of the input file (excluding ID).
## Example: if files are called /home/user/cohort_chr<1-22>.bim|bed|fam,
## set basename="/home/user/cohort_chr"
basename=$2

## Read path to the GENOSCORES analysis R script. Here you can use
## example.analysis.R
## Example: gscript="/home/user/example.analysis.R"
gscript=$3

echo "GENOSCORES analysis script: $gscript will be executed."

## Read chromosome number from command line
chr=$4

## Function to process one chromosome at a time.
## If you target file has suffixes after chromosome ID, for example:
## 'cohort_chr21_filtered', they need to be manually added to inputfile:
## inputfile="${basename}/cohort_chr${chr}_filtered"
processchrom() {
    inputfile="${basename}/RA_Control.final_${chr}"
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

processchrom

echo "All analyses have been executed."
