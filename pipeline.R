## This is the main script of the trans-QTL pipeline.
#' STEPS:
#' 1. Compute scores in the target cohort for eQTLs, pQTLs, other traits using
#'    GENOSCORES (see genoscores subdirectory).
#' 2. Aggregate regional scores into genome-wide trans scores
#'    (see gw.trans.scores.R).
#' 3. Load computed scores (see qtls.metadata.R)
#' 4. Load data on binary or continuous trait of interest.
#'    For example: T1D, or age of onset of T1D. See summary.stats.R
#' 5. Compute associations between the scores and the phenotype of interest.
rm(list=ls())

library(data.table)
library(GenomicRanges)
library(readxl)
library(ggplot2)
library(extraDistr)
library(genoscores)
library(ggcorrplot)
library(ggrepel)

## do we include HLA region in association analysis? We set it to FALSE in T1D
include.hla <- FALSE

## association FDR threshold
fdr=0.01

flanksize.kb <- 200 # used to create range for top SNP association in this dataset
DMhits.flanksize.kb <- 200 ## used to create GRanges object for reported T1D and T2D SNPs
regions.maxsep <- 1E5 ## max 100 kb between regions

## set to TRUE if you want to compute associations from scratch
newrun=TRUE

## STEP 1: Compute scores in the target cohort (genoscores subdirectory).
#'
#' - Configure `example.analysis.R` manually (details in genoscores/readme.md):
#'   - set `studyids` for which scores are to be computed.
#'     see https://genoscores.cphs.mvm.ed.ac.uk/studies for the list of
#'     available studies. For example: `studyid<-36` for eQTLGen.
#'   - check that the ancestry is EUR, or set a different ancestry if required.
#'   - set `ld.refplinkfile` path.
#' - Specify input files and output directories in `analysis.cli.sh`
#'   - check that `inputfile` in `analysis.cli.sh` respects your target
#'     file naming convention.
#'   @param outputdir File path to the directory where GENOSCORES will store
#'          computed scores.
#'   @param basename Base name (with full path) of the input file (excluding
#'          ID).
#'   @gscript Set path to the GENOSCORES analysis R script (example.analysis.R).
#'   @note You can create multiple example scripts - each for a different study,
#'         or you can provide a vector of many studies to `example.analysis.R`.
outputdir <- file.path(..., "eqtls")
basename <- file.path(...)
gscript <- file.path("genoscores", "example.analysis.R")

cmd <- sprintf("bash analysis.cli.sh %s %s %s", outputdir, basename, gscript)
system(cmd)

## STEP 2: Aggregate scores into cis, trans abd genome-wide trans scores.
#'
#' @param analysis String specifying the type of analysis. We currently support
#'        `eQTL`, `pQTL` (limited to on server computation), `immunecells`.
#'        other types of analysis can also be accommodated by modifying
#'        `gw.trans.scores.R` script.
#' @param base.in.dir Full path to the directory with GENOSCORES results. This
#'        should typically be the same as `outputdir` above.
#' @param base.out.dir Full path to the directory where the aggregated scores
#'        will be saved.
analysis <- ...
base.in.dir <- outputdir
base.out.dir <- file.path(...)

source("gw.trans.scores.R")

## STEP 3: Load computed scores and process metadata.
#'
#' @param base.out.dir Full path to the directory where the aggregated scores
#'        will be saved.
#' @source helperfunctions.R is a script with helper functions used to perform
#'         low-level operations.
#' @source diabetesgenes.R loads a list of known diabetes genes which will be
#'         used for annotation.
#' @source qtls.metadata.R contains functions `load.eqtl.data` and
#'         `load.pqtl.data` which are used to load and process eQTL and pQTL
#'         scores
source("helperfunctions.R")
source("diabetesgenes.R")
source("qtls.metadata.R")

load.eqtl.data(base.out.dir)
# load.pqtl.data(base.out.dir, genes.info)

pheno.file <- file.path(...)

## TODO: write guidelines in the script
source("phenotype.R")

## TODO: clean the code and check that scores are annotated correctly
source("eqtls.R")
source("pqtls.R")

## TODO: add further steps and an example report

## show memory use
objmem <- 1E-6 * sapply(ls(), function(x) {object.size(get(x))})
objclass <- sapply(ls(), function(x) {class(get(x))})
obj.dt <- data.table(object=objects(), memory=objmem, class=objclass)
setorder(obj.dt, -memory)
obj.dt[1:5, ]
gc()


