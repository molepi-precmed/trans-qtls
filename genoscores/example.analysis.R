#' This is a template script with example calls to GENOSCORES package
#' functions. It can be used as a guide to starting your own analysis script.
#' Here we provide a set of steps in a correct order to compute genotypic
#' scores using genotypic data provided by the user along with the weights
#' from several published GWAS stored in the GENOSCORES database.
#'
#' This script can be sourced directly into R session or run using command
#' line interface (CLI) via Rscript supplying command line arguments.
#' The option to run this script via CLI may be useful to automate score
#' computation if the target cohort is split by chromosome. Refer to
#' analysis.cli.sh script at https://genoscores.cphs.mvm.ed.ac.uk/downloads for
#' details on how to use CLI with GENOSCORES.
##===========================================================================

##---------------------------------------------------------------------------
#' STEP 1: Registration and obtaining access permission
##---------------------------------------------------------------------------
#' This package uses an API to make web requests to the GENOSCORES database
#' of published GWAS and retrieve the corresponding weights.
#'
#' To use this functionality, please create a user account at
#' https://genoscores.cphs.mvm.ed.ac.uk/register
#'
#' You will need to provide a valid email address and create a user password.
#' Please do not use an existing password but create a new one.
#' Please store your user credentials safe and accessible as you will need
#' them to setup the GENOSCORES package.
#'
#' After the registration is complete, you will be able to setup the
#' GENOSCORES package locally, which is explained in STEP 2.

##---------------------------------------------------------------------------
#' STEP 2: Initialize GENOSCORES on your working computer
##---------------------------------------------------------------------------
#' To use GENOSCORES you need to authorize the package to make web requests
#' to the GENOSCORES web application on your behalf.
#'
#' To start the initialisation process (assuming the package has been
#' installed) in an R session call initialise.genoscores(). This will start an
#' interactive initialisation session. Follow the instructions.
#'
#' Upon successful initialisation, a hidden file .genoscores.yml will be
#' created in your Linux HOME directory. It contains some sensitive
#' information about your access to the online database. Please do not move,
#' copy or otherwise expose this information to any third parties.

##---------------------------------------------------------------------------
#' STEP 3: Specify parameters of the analysis
##---------------------------------------------------------------------------
#' 1. Set the output directory (the directory where the analysis data files
#'    will be saved).
#' 2. Specify studies and/or trait types which are to be used in the analysis.
#'    GENOSCORES database contains an extensive collection of published GWAS
#'    which can be used in computing regional scores. You can refer to the
#'    current list of studies in the database at
#'    https://genoscores.cphs.mvm.ed.ac.uk/studies
#' 3. Configure system resources - set the number of scores and register
#'    parallelisation, set the memory allocation for PLINK.
#' 4. Set parameters for input data processing.
#' 5. Set parameters for score computation and adjustment for linkage
#'    disequilibrium (LD).
#'

##---------------------------------------------------------------------------
#' STEP 4: Run the analysis
##---------------------------------------------------------------------------
#' GENOSCORES performs the following steps:
#'
#' 1. Extracts weights for chosen studies and/or trait types from the
#' database corresponding to the SNPs in the plink file, and prepares a
#' catalog object.
#'
#' Catalog object has a weight for each SNP consisting of:
#'  * name of the trait reported in the study
#'  * trait type (see https://genoscores.cphs.mvm.ed.ac.uk/studies for trait types)
#'  * tissue (as measured in the study)
#'  * chromosome
#'  * genomic position
#'  * effect allele
#'  * linear regression slope (beta or actual weight)
#'  * p-value
#'  * alternative allele
#'  * pubmedid of the study from which the current weight is taken
#'  * number of samples used in the study
#'  * number of cases and controls (if case-control study)
#'  * ancestry of the study
#'
#' The Catalog object is then saved as catalog.Rdata.gz in the output.dir
#'
#' 2. Uses the weights to compute genotypic scores which are returned as a
#' matrix and saved as genotypicscore.Rdata.gz file in the output.dir
#'
#' 3. Information on the scores is saved to scoresinfo.Rdata.gz file, and
#' are also written out to file scoresinfo.csv which is
#' used in downstream analyses (e.g. score-score correlation matrices).
#'
## --------------------------------------------------------------------------

## Read command line arguments and choose the mode of analysis execution
#' @note When running via CLI, the script expects two command line arguments
#'       to be supplied to Rscript in this order: output.dir inputfile
args <- commandArgs(trailingOnly=TRUE)
mode <- ifelse(length(args) == 2, "CLI mode", "R session mode")
cat("Running analysis in", mode, "\n")

## Load GENOSCORES package
library("genoscores")

##---------------------------------------------------------------------------
#' 3.1 Set output directory
##---------------------------------------------------------------------------
#' @param output.dir Directory where GENOSCORES will write the results. If
#'        output.dir is NULL, the current directory will be used.
#'
output.dir <- NULL ## is passed from the calling script


##---------------------------------------------------------------------------
#' 3.2 Specify which studies and/or trait types to use
##---------------------------------------------------------------------------
#' You can select which GWAS to be used in the analysis by providing a
#' vector of study ids and/or a vector of trait types. As a study may have
#' been conducted on different ancestries, the ancestries of interest must be
#' always specified through their ancestry code (such as "EUR").
#'
#' The complete list of studies, ancestry codes and trait types is available at
#' https://genoscores.cphs.mvm.ed.ac.uk/studies
#'
#' @param ancestries Sets the ancestry codes for which gwaids are to be returned
#'        (compulsory argument)
#' @param studyids Sets which studies are to be used in the analysis
#' @param trait.types Sets which types of traits are to be used
#'
#' @note Here we use the following studies as an example:
#'          * Rheumatoid arthritis meta-analysis (Okada 2014)
#'          * Inflammatory bowel disease (Liu 2015)
#'          * Type 1 diabetes meta-analysis (Barrett 2009)
#'
studyids <- c(1, 2, 4)
ancestries <- "EUR"
trait.types <- NULL

##---------------------------------------------------------------------------
#' 3.3 Configure system resources
##---------------------------------------------------------------------------
#' @param num.cores Sets how many cores GENOSCORES can use on user's computer
#' @param num.batches Multiplication of the genotype matrix by the matrix of
#'        weights is performed on the user's computer. It is typically a
#'        resource-consuming operation and can be done in batches to
#'        reduce memory consumption. Default number of batches is 10, but it
#'        can be chosen per the available computing resources.
#' @param memory.plink Amount of memory reserved by plink in megabytes. If
#'        \code{NULL}, plink tries to reserve half of the system's RAM for
#'        its workspace.
#'
num.cores <- 5
num.batches <- 15
memory.plink <- 20000

## Register parallelisation mode if more than one core is provided
options(cores=num.cores)

##---------------------------------------------------------------------------
#' 3.4 Set parameters for input data processing
##---------------------------------------------------------------------------
#' @param inputformat Set input file format.
#'        The input file should contain genotype data from the patients for
#'        which the regional genotypic scores are to be computed.
#'        We currently support "plink", "dosage", "bgendosage" as input file
#'        formats. Note, that we currently support the input file
#'        pre-processing only for "plink" files.
#' @param inputfile Specify the input file location.
#' @param prepare.plinkfile Prepare input PLINK file by updating SNP ids in
#'        the input file to the most recent build.
#' @param filter.snpids Filter SNP ids in the input file to keep only those
#'        present in the GENOSCORES database.
#' @param liftover.file liftOver chain file used to lift the genome positions
#'        to build hg38. You can download various chain files
#'        from https://hgdownload.soe.ucsc.edu/downloads.html
#'        If your build is already hg38 set this variable to NULL. To find
#'        out the build of your input file call guess.genome.build(input.file)
#'        first.
#'
inputformat <- "plink"
inputfile <- NULL ## will be passed from the calling script
prepare.plinkfile <- TRUE
filter.plinkfile <- TRUE
liftover.file <- NULL


##--------------------------------------------------------------------------
#' 3.5 Set parameters for score computation and adjustment for LD
##--------------------------------------------------------------------------
#' @param gwasids.region Allows to control which \code{gwasids} are split into
#'        regional scores, splitting SNPs from each GWAS into regions. It can
#'        be set to \code{"all"} (default) or \code{"none"}, or to a numeric
#'        vector that specifies the \code{gwasids} to be considered complex.
#' @param ld.refplinkfile Path to a reference plink file to be used for
#'        calculating the LD and adjusting the weights (GWAS betas) for
#'        LD structure. If NULL, no LD adjustment is performed.
#'        It is strongly recommended to provide a file with the genotypes of a
#'        reference population (such as 1000G or UK10K) which matches the
#'        ancestry of the cohort being analysed. For more details please refer
#'        to the online documentation and the documentation of get.scores().
#' @param filter.pvalue Set a threshold on p-value to filter SNPs.
#'        All SNPs at p < filter.pvalue will be discarded from the analysis.
#' @param region.pvalue Genome is split into regions for which regional
#'        scores are computed. A region is created when there is at least one
#'        SNP with p-value that satisfies the criterion p < region.pvalue.
#' @param region.distance Distance (in base pairs) above which two subsequent
#'        SNPs are split into different regions.
#' @param ld.method The type of LD adjustment to perform.
#'        Available methods: "pseudoinverse" (default), "ppca", "ppca.adj",
#'        "addlambda.diag".
#' @param ld.eigvar.thresh A numeric scalar in [0, 1] specifying the
#'        percentage of variance we want to retain when a pseudoinverse is
#'        used to invert the LD matrix.
#' @param ld.prune.thresh A scalar threshold in [0, 1] for pruning SNPs
#'        in very high LD, or NULL to avoid pruning.
#'
gwasids.region <- "all"
ld.refplinkfile <- file.path(...) # or NULL for no LD correction (discouraged)
opts <- genoscores.opts(filter.pvalue=1E-5,
                        region.pvalue=1E-5,
                        region.distance=5E5,
                        ld.method="pseudoinverse",
                        ld.eigvar.thresh=0.9,
                        ld.prune.thresh=NULL)

##--------------------------------------------------------------------------
#' 4. Notes on running the analysis
##--------------------------------------------------------------------------
#' Computation of genotypic scores for any real-life genotypic datasets is
#' a resource-intensive task. Therefore, we recommend running GENOSCORES on
#' a multiple cores machines with sufficient memory.
#'
#' As an example the analysis with the full genome data for ~500 individuals
#' requesting weights from 3 autoimmune diseases on 4 cores machine with
#' 28GB of memory took around 20 minutes of user time.
#'
#' It is also advised to limit the memory consumed by the GENOSCORES to
#' prevent system crashed. This can be achieved via the following command:
#' ulimit -v 30000000
#' Here, we allow the script to consume no more than 30GB of memory. It is
#' recommended to leave at least few GB of memory and at least one code to
#' the correct operation of the Linux OS while running the analysis.
#'

#################   usually no need to edit below this line   ###############

## settings extracted from the command line
if (mode == "CLI mode") {
  output.dir <- args[1]
  inputfile <- args[2]
}



cat("Downloading gwasids for requested studies from the database ... \n")
gwasids <- select.gwasids(ancestries=ancestries, studyids=studyids,
                          trait.types=trait.types)

## prepare PLINK file
if (inputformat == "plink") {
    if (prepare.plinkfile) {
        cat("Preparing plink file ...\n")
        ptm <- proc.time()
        prepare.plinkfile(plinkfile=inputfile,
                          liftover.file=liftover.file,
                          output.dir=output.dir,
                          memory.plink=memory.plink)
        ptm <- proc.time() - ptm
        cat("Finished in", ptm[3], "seconds.\n")
    }
} else if (prepare.plinkfile) {
    cat("Update rsid functionality not implemented for", inputformat, "format.")
}

## filter PLINK file
if(inputformat == "plink") {
    if (filter.plinkfile) {
        cat("Filtering rsids ...\n")
        ## If "prepare PLINK file" step was run, then we want to use the
        ## file outputed by PLINK in place of the user specified inputfile
        inputfile.updated <- file.path(output.dir, paste0(basename(inputfile),
                                                          "_snpids.updated"))
        inputfile <- ifelse(file.exists(paste0(inputfile.updated, ".bim")),
                            inputfile.updated, inputfile)
        ptm <- proc.time()
        filter.plinkfile(plinkfile=inputfile,
                         memory.plink=memory.plink,
                         output.dir=output.dir)
        ptm <- proc.time() - ptm
        cat("Finished in", ptm[3], "seconds.\n")
    }
} else if (filter.plinkfile) {
    cat("Filter rsid functionality not implemented for", inputformat, "format.")
}

## Choose which file should we use as input to compute scores.
## Since it is recommended to compute scores using file which has been
## correctly preprocessed, we print a warning if such file cannot be found.
files <- list.files(output.dir, pattern="_snpids.updated|_forscores",
                    full.names=TRUE)
regex <- c("_snpids.updated_forscores", "_forscores", "_snpids.updated")
matchedfile <- unlist(sapply(regex, function(z) grep(z, files, value=TRUE)))[1]
inputfile <- ifelse(is.na(matchedfile), inputfile,
                    gsub(".bed$|.bim$|.fam$", "", matchedfile))
cat("Selected", inputfile, "for the analysis\n")

## compute scores
cat("Calculating genotypic scores ...\n")
ptm <- proc.time()
scores <- get.scores(inputfile=inputfile,
                     format=inputformat,
                     gwasids=gwasids,
                     gwasids.region=gwasids.region,
                     ld.refplinkfile=ld.refplinkfile,
                     opts=opts,
                     output.dir=output.dir,
                     num.batches=num.batches)
ptm <- proc.time() - ptm
cat("Finished score computation with", num.cores, "cores in:\n")
print(ptm)
