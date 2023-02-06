#' Write code to load and process your phenotype here. Below code is a basic
#' example of how to read the phenotype file and ensure exact match between the
#' samples samples in the scores and in the phenotype.
## =============================================================================
## read individual-level phenotype dataset
phenotype <- fread(pheno.file)

## by default GENOSCORES merges the family id (FID) and individual id (IID)
## into a single identifier: FID_IID. You may need to ensure that sample ids in
## the phenotype file respect this convention. Refer to PLINK .fam file specs
## for more details on IDs: https://www.cog-genomics.org/plink2/formats#fam
## For example:
phenotype[, id := paste0(ID_1, "_", ID_2)]

## you can check that the sample ids are matched (in correct order) between the
## phenotype and the scores
phenotype <- phenotype[match(rownames(genome.wide.scores), id), ]

## ensure no missing data in phenotypes
phenotype <- phenotype[complete.cases(phenotype)]

## if some samples were removed from samples, we also remove them from scores
sampleids <- intersect(phenotype[, id], rownames(genome.wide.scores))
genome.wide.scores <- genome.wide.scores[rownames(genome.wide.scores) %in% sampleids, ]

## ensure sample ids are identical
stopifnot(identical(phenotype[, id], rownames(genome.wide.scores)))
