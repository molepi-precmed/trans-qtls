## CCL15 represented by 2 cis scores in pqtl.cscores.info X_632153_1_17 and X_632154_1_17
## for two different GWAS ids.
## but only one column labelled CCL15 in pqtl.scores
## probably this is the sum over 2 GWAS studies
pqtl.cscores.info[gene_symbol=="CCL15"]
grep("CCL15", colnames(pqtl.scores))

pqtl.scoreids <- colnames(pqtl.scores)

pqtl.scores.dt <- data.table(id=rownames(pqtl.scores), pqtl.scores)
setkey(pqtl.scores.dt, id)
cpep.pqtl.scores <- merge(pqtl.scores.dt, cpep.dt, all=TRUE)
aod.pqtl.scores <- merge(pqtl.scores.dt, aod.dt, all=TRUE)

## the above merge for scores in t1bio cohort will result in NAs for the
## controls which we need to filter out at this stage
cpep.pqtl.scores <- cpep.pqtl.scores[complete.cases(cpep.pqtl.scores),] ## 4958
aod.pqtl.scores <- aod.pqtl.scores[complete.cases(aod.pqtl.scores)] ## 5342

## scale continuous phenotype
cpep.pqtl.scores[, duration.cpeptide := scale(duration.cpeptide)]
aod.pqtl.scores[, PHENO := scale(PHENO)]

## calculate SD of each score and merge with metadata.
## this is not strictly correct as we are merging the SD of the genomewide trans score with the first locus-specific score for that gene symbol
pqtl.scores.sd <- apply(pqtl.scores.dt[, -1], 2, sd)
pqtl.scores.sd <- data.table(matrix.colname=names(pqtl.scores.sd),
                             sdscore=as.numeric(pqtl.scores.sd))

## FIXME: this creates a new column named gene_symbol.1
setkey(pqtl.scores.sd, matrix.colname)
setkey(pqtl.allscores.info, matrix.colname)
pqtl.allscores.info <- pqtl.scores.sd[pqtl.allscores.info]

# pqtls.file <- "pqtl.coeffs.RData"
pqtls.file <- "cpep.pqtl.coeffs.RData"
pqtls.file <- "aod.pqtl.coeffs.RData"
if(newrun) {
    pqtl.coeffs <- NULL
    for(scoreid in pqtl.scoreids) {
        score.model <- glm(formula=PHENO ~ get(scoreid) + PC1 + PC2 + PC3,
                           family="gaussian",
                           data=aod.pqtl.scores)
        coeff <- summary(score.model)$coefficients[2, , drop=FALSE]
        rownames(coeff) <- scoreid
        pqtl.coeffs <- rbind(pqtl.coeffs, coeff)
    }
    save(pqtl.coeffs, file=pqtls.file)
} else {
    load(pqtls.file)
}

pqtl.coeffs <- data.table(matrix.colname=rownames(pqtl.coeffs), pqtl.coeffs)
colnames(pqtl.coeffs)[5] <- "pvalue"

setkey(pqtl.coeffs, matrix.colname)
setkey(pqtl.allscores.info, matrix.colname)
pqtl.coeffs <- pqtl.allscores.info[pqtl.coeffs]

#pqtl.hits <- pqtl.coeffs[minuslog10p > 6]
#pqtl.hits[, .(matrix.colname, gene_symbol, numgwas, locus.diversity,
#              gene_chrom, gene_startpos,
#              qtl_type, numscores, Estimate, minuslog10p)]

pqtl.coeffs[, cisx := ifelse(qtl_type=="cis-x", 1, 0)]
pqtl.coeffs[qtl_type=="cis-x", qtl_type := "cis"]

pqtl.coeffs[, minuslog10p := -log10(pvalue)]

setnames(cpep.pqtl.coeffs, old="Std. Error", new="SE")
setnames(aod.pqtl.coeffs, old="Std. Error", new="SE")

setorder(pqtl.coeffs, -minuslog10p)
pqtl.coeffs.unique <- unique(pqtl.coeffs, by=c("gene_symbol", "qtl_type"))

## cast pqtl.coeffs to wide format -- one row per gene_symbol, cis- and trans- effects
pqtl.coeffs.wide <- dcast(pqtl.coeffs.unique,
                     gene_symbol + gene_chrom + gene_startpos + gene_endpos ~ qtl_type,
                     value.var=c("numscores", "Estimate", "pvalue", "locus.diversity",
                                 "minuslog10p"))
pqtl.coeffs.wide[, gene_chrom := factor(as.character(gene_chrom), levels=c(1:22, "X", "Y"))]
setnames(pqtl.coeffs.wide, "gene_chrom", "chromosome")

setkey(genes.info, gene_symbol)
setkey(pqtl.coeffs.wide, gene_symbol)

## merge cpep.coeffs on gene_symbol with genes.info, which contains boolean columns for is.t1dgene and is.t2dgene
pqtl.coeffs.wide <- genes.info[, .(gene_symbol, gene_biotype, is.t1dgene, is.t2dgene,
                            near.t1dgenes, near.t2dgenes)][pqtl.coeffs.wide]


pqtl.topscores.colnames <- pqtl.coeffs[minuslog10p > 6, matrix.colname]
