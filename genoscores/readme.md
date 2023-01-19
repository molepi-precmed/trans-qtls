# Compute regional scores in a target cohort

We compute the regional scores in the specified target cohort using weights from GENOSCORES database.

1. [example.analysis.R](example.analysis.R) is an example R script containing the instructions for computing scores. Please reas carefully to understand the procedure. The things which have to be modified in this script are:
  - `studyids` -- see https://genoscores.cphs.mvm.ed.ac.uk/studies for the list of available studies. For example: `studyid<-36` for eQTLGen.
  - `ancestries` set to `EUR` unless using a multi-ncestry study. 
  - `ld.refplinkfile` specify a full path to the LD reference file in PLINK format. If required, we can provide the latest release of 1000 Genomes in hg38.
2. [analysis.cli.sh](analysis.cli.sh) is a BASH script which executes  [example.analysis.R](example.analysis.R) batching analysis per chromosome. Here, you need to specify the following:
  - `outputdir` -- directory to store scores (once subdirectory per chromosome will be created automatically).
  - `basename` -- name (with full path) of the input file (see script for explanation).
  - `gscript` -- path to the GENOSCORES script.
  - `N` -- how many chromosomes do you want to run in parallel? This depends on the size of the target cohort and available computational resources.
  - `inputfile` -- modify accordingly to the full name of your target cohort files. If the names have any suffixes after chromosome ID, these need to be added.
  
## General notes

1. If multiple analyses are to be run (eQTL, pQTL, etc), then it may be helpful to either:
   - specify multiple study ids in [example.analysis.R](example.analysis.R)
   - write a wrapper script around [analysis.cli.sh](analysis.cli.sh) which will call different versions of [example.analysis.R](example.analysis.R), each specifying its own study. This approach may be beneficial if, say, eQTL and pQTL analyses require different thresholds.
