## Main script of the trans-QTL pipeline
##==============================================================================
#' This is a template which sources the scripts responsible for the key steps
#' of computing regional genotypic scores, aggregating the genome wide trans
#' scores and computing associations between these scores and the phenotype of
#' interest.
#' This template can be easily modifiable to accomodate the requirements of your
#' analysis. For details, see each sourced script and the corresponding READMEs
#' or comments in code.
#'
#' Key steps:
#' 1. Compute scores in the target cohort for selected intermediate traits.
#'    Example: blood eQTLs, plasma pQTLs, etc.
#'    See genoscores subdirectory
#' 2. Aggregate regional scores into genome-wide trans scores.
#' 3. Load computed scores.
#' 4. Load data on binary or continuous trait of interest.
#'    For example: T1D, or age of onset of T1D.
#' 5. Compute associations between the scores and the phenotype of interest.
##==============================================================================
rm(list=ls())

library(data.table)
library(GenomicRanges)
library(readxl)
library(ggplot2)
library(extraDistr)
library(genoscores)
library(ggcorrplot)
library(ggrepel)
library(parallel)
library(foreach)

source("helperfunctions.R")

##------------------------------------------------------------------------------
#' STEP 0: Specify generic parameters
##------------------------------------------------------------------------------
#' @param include.hla Logical specifying if HLA region is included in the
#'                    analysis. We set it to FALSE in T1D analysis.
#' @param fdr Double specifying FDR threshold for associations.
#' @param flanksize.kb Integer (currently not used) create range for top SNP
#'                     association in this dataset.
#' @param DMhits.flanksize.kb Integer create GRanges object for reported T1D
#'                            and T2D SNPs.
#' @param regions.maxsep Integer (currently not used) max 100 kb between
#'                       regions.
#' @param newrun Logical setting whether to compute associations from scratch.

include.hla <- FALSE
fdr=0.01
flanksize.kb <- 200
DMhits.flanksize.kb <- 200
regions.maxsep <- 1E5

newrun=TRUE

##------------------------------------------------------------------------------
#' STEP 1: Compute scores in the target cohort
##------------------------------------------------------------------------------
#' Uses GENOSCORES package to query GWAS weights for the specified study and
#' compute regional scores in the target cohort, which are then adjusted for LD
#' using specified LD reference cohort.
#' Check https://genoscores.cphs.mvm.ed.ac.uk/quickstart to learn how to
#' install and use the package.
#'
#' The example/helper scripts to get you started with computing the scores as
#' well as a more detailed README are located in `genoscores` subdirectory:
#' - `genoscores/example.analysis.R` is an example R script containing the
#'    instructions for computing scores. Please read carefully to understand
#'    the procedure.
#' - is a BASH script which executes `example.analysis.R` batching analysis per
#'   chromosome.
#'
#' The minimal changes to `genoscores/example.analysis.R` and
#' `genoscores/`analysis.cli.sh` scripts required to compute the scores are:
#' 1. set `studyids` in `genoscores/example.analysis.R`. These specify studies
#'    for which scores are to be computed. Refer to
#'    https://genoscores.cphs.mvm.ed.ac.uk/studies for the list of studies.
#'    Example: `studyid<-36` for eQTLGen.
#' 2. check that the ancestry is `EUR` in `genoscores/example.analysis.R`,
#'    or set a different ancestry if required.
#' 3. set `ld.refplinkfile` path in `genoscores/example.analysis.R`. This should
#'    be the path to LD reference PLINK file in hg38. We can provide the most
#'    recent release of 1000 Genomes if needed.
#' 4. Specify input file and output directory in `genoscores/analysis.cli.sh`.
#'    Check that `inputfile` in `analysis.cli.sh` respect your target file
#'    naming convention.
#'
#' Example of running `analysis.cli.sh` is included below.
#' @param output.dir File path to the directory where GENOSCORES will store
#'        computed scores.
#' @param basename Base name (with full path) of the input file (excluding ID).
#' @param gscript Set path to the GENOSCORES analysis R script.
#' @note You can create multiple example scripts - each for a different study,
#'         or you can provide a vector of many studies to `example.analysis.R`.

output.dir <- file.path(...)
basename <- file.path(...)
gscript <- file.path("genoscores", "example.analysis.R")

cmd <- sprintf("bash genoscores/analysis.cli.sh %s %s %s",
               output.dir, basename, gscript)
system(cmd)

##------------------------------------------------------------------------------
#' STEP 2: Load regional scores, classify into cis, cis-x and trans and
#' aggregate trans scores into genome-wide trans scores.
##------------------------------------------------------------------------------
#' @param analysis String specifying the type of analysis. We currently support:
#'        - "eQTL"
#'        - "pQTL" (limited to on server computation)
#'        - "immunecells"
#'        you can implement other traits by modifying `gw.trans.scores.R`.
#' @param score.dir Full path to the directory with GENOSCORES results. This
#'        should typically be the same as `outputdir` above.
#' @param output.dir Full path to the directory where the aggregated scores
#'        will be saved.
analysis <- ...
score.dir <- output.dir
output.dir <- file.path(...)

source("gw.trans.scores.R")

##------------------------------------------------------------------------------
#' STEP 3: Load computed genome-wide trans scores and the metadata. Optionally
#' create a table of annotations for genes implicated in the disease of interest
#' (T1D and T2D in this case).
##------------------------------------------------------------------------------
## (optional) load the script to create a table of annotated genes implicated in
## the disease of interest
source("diabetesgenes.R")

## load eQTL or pQTL scores and metadata depending on `analysis` variable
score.dir <- output.dir
source("qtl.metadata.R")

##------------------------------------------------------------------------------
#' STEP 4: Load binary of continuous phenotype and match samples between scores
#' and the phenotype.
##------------------------------------------------------------------------------
pheno.file <- file.path(...)

## script below is a short guideline. You may need to extend it according to the
## specific analysis requirements.
source("phenotype.R")

##------------------------------------------------------------------------------
#' STEP 5: Compute associations between the scores and the phenotype of
#' interest.
##------------------------------------------------------------------------------
#' @param binary Logical specifying if case/control (TRUE) or continuous (FALSE)
#'        phenotype is to be analysed.
#' @param output.dir Full path to the directory where the association results
#'        will be saved.
binary <- TRUE
output.dir <- file.path(...)

## register parallelisation mode to run several chromosomes at a time
options(cores=4)

source("association.R")

## show memory use
objmem <- 1E-6 * sapply(ls(), function(x) {object.size(get(x))})
objclass <- sapply(ls(), function(x) {class(get(x))})
obj.dt <- data.table(object=objects(), memory=objmem, class=objclass)
setorder(obj.dt, -memory)
obj.dt[1:5, ]
gc()
