#' Merge scores with the phenotype and test for association.
## =============================================================================
## convert scores matrix to data table
scoreids <- colnames(genome.wide.scores)
scores.dt <- data.table(id=rownames(genome.wide.scores), genome.wide.scores)

## stop with an error if scores and phenotype have different samples
stopifnot(identical(phenotype[, id], rownames(genome.wide.scores)))

## scale continuous phenotype
if (!binary) phenotype[, PHENO := scale(PHENO)]

## merge scores with phenotype
setkey(scores.dt, id)
phenotype.scores <- merge(scores.dt, phenotype, all=TRUE)

## code below runs a simple parallelised for loop to test each score for
## association with the phenotype according to the specified formula.
## We assumed that the column with phenotype values is called "PHENO".
## You can also add covariates similarly to how we added the principal
## components PC1, PC2 and PC3 below.
coeffs <- NULL
model.family <- ifelse(binary, "binomial", "gaussian")
if(newrun) {
    coeffs <- foreach(scoreid=scoreids, .combine=rbind) %dopar% {
        score.model <- glm(formula=PHENO ~ get(scoreid) + PC1 + PC2 + PC3,
                           family=model.family,
                           data=phenotype.scores)
        coeff <- summary(score.model)$coefficients[2, , drop=FALSE]
        rownames(coeff) <- scoreid

        coeff <- data.table(matrix.colname=rownames(coeff), coeff)
        colnames(coeff)[5] <- "pvalue"
        colnames(coeff)[3] <- "SE"
        if (binary) {
            coeff[, z.abs := abs(`z value`)]
            coeff[, minuslog10p := minus.log10.phnorm.extreme(z.abs)]
        } else {
            coeff[, minuslog10p := -log10(pvalue)]
        }
    }
    save(coeffs, file=file.path(output.dir, "coeffs.RData"))
} else {
    load(file.path(output.dir, "coeffs.RData"))
}

# Make sure that regressions were run for all scores
if(nrow(coeffs) != length(scoreids)){
    stop('Number of associations does not equal the number of scores.')
}

## merge coeffs with scores metadata
if (analysis=="eQTL")
    scoresinfo <- allscores.info
if (analysis=="pQTL")
    scoresinfo <- pqtl.allscores.info

setkey(coeffs, matrix.colname)
setkey(scoresinfo, matrix.colname)

cat(coeffs[!(matrix.colname %in% scoresinfo$matrix.colname), .N],
    "coeffs not matched in allscores.info\n")

coeffs <- scoresinfo[coeffs]
table(coeffs$qtl_type, exclude=NULL)

## All the post-processing steps below are optional. Some steps may be useful
## for converting the results to format which is easy to use in publication.

## manage cis, cis-x and trans labels
coeffs[is.na(qtl_type), qtl_type := ifelse(grep("trans$", matrix.colname), "trans", NA)]
coeffs[, cisx := ifelse(qtl_type=="cis-x", 1, 0)]
coeffs[qtl_type=="cis-x", qtl_type := "cis"]

## temporary fix for unmatched gene_symbol -- but these rows have regions NA
# coeffs[is.na(gene_symbol) & substr(matrix.colname, 1, 2) != "X_", gene_symbol := matrix.colname]

## compute standardised log odds ratio
if (binary) coeffs[, Estimate := Estimate * sdscore]

genes.duplicated <- coeffs[, .N, by=.(gene_symbol, qtl_type, cisx)][N > 1, gene_symbol]
cat(coeffs[gene_symbol %in% genes.duplicated, .N], "duplicated gene symbols\n")

## loci field is missing for cis-eQTLs
coeffs[is.na(regions), regions := "0"]

setorder(coeffs, -minuslog10p)
coeffs.unique <- unique(coeffs, by=c("gene_symbol", "qtl_type"))

## cast coeffs to wide format -- one row per gene_symbol, cis- and trans-
## effects
if (analysis == "eQTL") {
    coeffs.wide <- dcast(coeffs.unique,
                         gene_symbol + gene_chrom + gene_startpos + gene_endpos ~ qtl_type,
                         value.var=c("regions",
                                     "locus.diversity", "numscores", "Estimate", "pvalue",
                                     "minuslog10p"))
    setnames(coeffs.wide, "regions_trans", "regions")
    setnames(coeffs.wide, "locus.diversity_trans", "locus.diversity")
}
if (analysis == "pQTL") {
    coeffs.wide <- dcast(coeffs.unique,
                         gene_symbol + gene_chrom + gene_startpos + gene_endpos ~ qtl_type,
                         value.var=c("numscores", "Estimate", "pvalue", "locus.diversity",
                                     "minuslog10p"))
    setnames(coeffs.wide, "locus.diversity_trans", "locus.diversity")
}

setkey(coeffs.wide, gene_symbol)
setkey(genes.info, gene_symbol)

## merge coeffs.wide with annotation of known T1D and T2D genes
coeffs.wide <- genes.info[, .(gene_symbol, gene_biotype, is.t1dgene, is.t2dgene,
                              near.t1dgenes, near.t2dgenes)][coeffs.wide]

## annotate genes by whether they are known cause for monogenic diabetes
coeffs.wide[, monogenic := ifelse(gene_symbol %in% t1d.monogenes, "+", ".")]

save(coeffs.wide, file=file.path(output.dir, "coeffs.wide.RData"))
