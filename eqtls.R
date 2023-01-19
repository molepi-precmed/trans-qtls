## merge scores with case-control dataset and test for association
scoreids <- colnames(genome.wide.scores)
scores.dt <- data.table(id=rownames(genome.wide.scores), genome.wide.scores)
## calculate SD of each score and merge with metadata table allscores.info


setkey(scores.dt, id)
cc.scores <- merge(scores.dt, cc.dt, all=TRUE)

if(newrun) {
    coeffs <- NULL
    for(scoreid in scoreids) {
        score.model <- glm(formula=PHENO ~ get(scoreid) + PC1 + PC2 + PC3,
                           family="binomial",
                           data=cc.scores)
        coeff <- summary(score.model)$coefficients[2, , drop=FALSE]
        rownames(coeff) <- scoreid
        coeffs <- rbind(coeffs, coeff)
    }

    caseonly.coeffs <- NULL
    for(scoreid in scoreids) {
        caseonly.model <- lm(formula=get(scoreid) ~ hla.score + PC1 + PC2 + PC3,
                         data=cc.scores[PHENO==1])
        caseonly.coeff <- summary(caseonly.model)$coefficients[2, , drop=FALSE]
        rownames(caseonly.coeff) <- scoreid
        caseonly.coeffs <- rbind(caseonly.coeffs, caseonly.coeff)
    }
    save(coeffs, caseonly.coeffs, file="coeffs.RData")
} else {
    load("coeffs.RData")
}

coeffs <- data.table(matrix.colname=rownames(coeffs), coeffs)
colnames(coeffs)[5] <- "pvalue"
coeffs[, z.abs := abs(`z value`)]
coeffs[, minuslog10p := minus.log10.phnorm.extreme(z.abs)]

caseonly.coeffs <- data.table(matrix.colname=rownames(caseonly.coeffs), caseonly.coeffs)
colnames(caseonly.coeffs)[5] <- "pvalue"
caseonly.coeffs[, z.abs := abs(`t value`)]
caseonly.coeffs[, minuslog10p := minus.log10.phnorm.extreme(z.abs)]

## merge coeffs with allscores.info
## this will also merge region numbers
## but we want the scoreids for each gene
setkey(coeffs, matrix.colname)
setkey(allscores.info, matrix.colname)
cat(coeffs[!(matrix.colname %in% allscores.info$matrix.colname), .N],
    "coeffs not matched in allscores.info\n")

setkey(caseonly.coeffs, matrix.colname)
cat(caseonly.coeffs[!(matrix.colname %in% allscores.info$matrix.colname), .N],
    "caseonly coeffs not matched in allscores.info\n")

caseonly.coeffs[!(matrix.colname %in% eqtl.info$gene_symbol) &
       !(matrix.colname %in% eqtl.info$scoreid), matrix.colname]

coeffs <- allscores.info[coeffs]
table(coeffs$qtl_type, exclude=NULL)
coeffs[is.na(qtl_type), qtl_type := ifelse(substr(matrix.colname, 1, 2) == "X_",
                                           "trans", "cis")]
coeffs[is.na(qtl_type) & qtl_type=="cis", gene_symbol := matrix.colname]

setorder(coeffs, -minuslog10p)
#coeffs <- unique(coeffs, by=c("gene_symbol", "qtl_type"))

caseonly.coeffs <- allscores.info[caseonly.coeffs]

## temporary fix for unmatched gene_symbol -- but these rows have regions NA
#coeffs[is.na(gene_symbol) & substr(matrix.colname, 1, 2) != "X_", gene_symbol := matrix.colname]

coeffs[, chromosome := factor(chromosome, levels=c(1:22, "X", "Y"))]
coeffs[, Estimate := Estimate * sdscore] ## standardize log odds ratio

caseonly.coeffs[, chromosome := factor(chromosome, levels=c(1:22, "X", "Y"))]
caseonly.coeffs[, Estimate := Estimate * sdscore] ## standardize log odds ratio

## QQ plot
p.cc.qq <- ggplot(data=coeffs, aes(sample=z.abs)) +
               ggplot2::stat_qq(distribution=extraDistr::qhnorm) +
        ggplot2::stat_qq_line(distribution=extraDistr::qhnorm) +
        xlab("Quantile of null (standard half-normal) distribution") +
        ylab("Quantile of absolute value of test statistic")
#p.cc.qq + scale_y_continuous(sec.axis=sec_axis(minus.log10.phnorm.extreme,
#                                               name="minus log10 p-value"))
#p.cc.qq

p.threshold.cc <- p.threshold(10^-coeffs$minuslog10p)

genes.duplicated <- coeffs[, .N, by=.(gene_symbol, qtl_type, cisx)][N > 1, gene_symbol]
cat(coeffs[gene_symbol %in% genes.duplicated, .N], "duplicated gene symbols\n")

## loci field is missing for cis-eQTLs
coeffs[is.na(regions), regions := "0"]


setorder(coeffs, -minuslog10p)
coeffs.unique <- unique(coeffs, by=c("gene_symbol", "qtl_type"))

## cast coeffs to wide format -- one row per gene_symbol, cis- and trans- effects
coeffs.wide <- dcast(coeffs.unique,
                     gene_symbol + chromosome + gene_startpos + gene_endpos ~ qtl_type,
                     value.var=c("regions", # "t1deQTL.nums",
                                 "locus.diversity", "numscores", "Estimate", "pvalue",
                                 "minuslog10p", "z.abs"))
setkey(coeffs.wide, gene_symbol)
setnames(coeffs.wide, "regions_trans", "regions")
setnames(coeffs.wide, "locus.diversity_trans", "locus.diversity")
#setnames(coeffs.wide, "t1deQTL.nums_trans", "t1deQTL.nums")

######################################

setorder(caseonly.coeffs, -minuslog10p)
caseonly.coeffs.unique <- unique(caseonly.coeffs, by=c("gene_symbol", "qtl_type"))

## cast caseonly.coeffs to wide format -- one row per gene_symbol, cis- and trans- effects
caseonly.coeffs.wide <- dcast(caseonly.coeffs.unique,
                     gene_symbol + chromosome + gene_startpos + gene_endpos ~ qtl_type,
                     value.var=c("regions", # "t1deQTL.nums",
                                 "locus.diversity", "numscores", "Estimate", "pvalue",
                                 "minuslog10p", "z.abs"))
setkey(caseonly.coeffs.wide, gene_symbol)
caseonly.coeffs.wide <- caseonly.coeffs.wide[, .(gene_symbol,
                                                 case.only.effect_cis=Estimate_cis,
                                                 case.only.minuslog10p_cis=minuslog10p_cis,
                                                 case.only.effect_trans=Estimate_trans,
                                                 case.only.z.abs_trans=z.abs_trans,
                                                 case.only.minuslog10p_trans=minuslog10p_trans)]
coeffs.wide <- caseonly.coeffs.wide[coeffs.wide]

## merge annnotations  --
## t1d.genes, t2d.genes in transcript region, neargenes in source regions --
## with coeffs.wide
setkey(genes.info, gene_symbol)
setkey(coeffs.wide, gene_symbol)

## merge cc.coeffs on gene_symbol with genes.info, which contains boolean columns for is.t1dgene and is.t2dgene
cc.coeffs <- genes.info[, .(gene_symbol, gene_biotype, is.t1dgene, is.t2dgene,
                            near.t1dgenes, near.t2dgenes)][coeffs.wide]

cc.coeffs[, gene_startpos := gene_startpos]
cc.coeffs[, pvalue_cis := format.z.aspvalue(z.abs_cis)]
cc.coeffs[, pvalue_trans := format.z.aspvalue(z.abs_trans)]
cc.coeffs[, case.only.pvalue_trans := format.z.aspvalue(case.only.z.abs_trans)]

cc.coeffs[, pvalue_cis := ifelse(minuslog10p_cis > 2, format.z.aspvalue(z.abs_cis), ".")]
cc.coeffs[, Estimate_cis := ifelse(minuslog10p_cis > 2, Estimate_cis, NA)]

cc.coeffs[, monogenic := ifelse(gene_symbol %in% t1d.monogenes, "+", ".")]

