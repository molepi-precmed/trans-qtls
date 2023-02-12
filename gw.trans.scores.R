#' Process genotypic scores computed for the target cohort using various traits
#' and aggregate them into genome-wide trans scores. This script assumes that
#' scores were computed using `genoscors/analysis.cli.sh` and were saved into
#' one subdirectory per chromosome.
#' This script expects `analysis`, `score.dir`, `output.dir`, `include.hla`
#' exist in the global environment.
#'
#' This is just an example script which can be run as is or modified/extended
#' to cover other types of scores (for example mQTLs).
## =============================================================================

#' Read in the scores matrix and metadata from each chromosome subdirectory
##------------------------------------------------------------------------------
all.scores <- NULL
all.scoresinfo <- NULL
for (idx in 1:22) {
    cat("Processing scores on chromosome", idx, "\n")
    load(file.path(score.dir, idx, "genotypicscore.Rdata.gz"))
    all.scores <- cbind(all.scores, score)
    cat("Processing scoresinfo on chromosome", idx, "\n")
    load(file.path(score.dir, idx, "scoresinfo.Rdata.gz"))
    all.scoresinfo <- rbind(all.scoresinfo, scoresinfo)
}
rm(score)
rm(scoresinfo)

#' Perform some post-processing of scores and metadata. This can be modified
#' according to the requirements of the specific analyses.
##------------------------------------------------------------------------------

## we usually remove background scores which were not assigned to any region and
## contain mostly noise
all.scoresinfo <- all.scoresinfo[region!=0]

## if HLA region is not included in the analysis, we remove any scores which
## intersect it.
## HLA region in HG38: chr6:28,510,120-33,480,577
## https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC
## here we use a slightly lower startpos to exclude any genes which may
## be stretching out from the HLA region.
if (include.hla) {
    ans <- all.scoresinfo
} else {
    hla <- data.table(startpos=25000000, endpos=34000000)
    setkey(hla, startpos, endpos)

    ## make an overlap join (join by range) between scoresinfo and HLA region
    idx <- foverlaps(all.scoresinfo, hla, type="within", which=TRUE)

    ## extract scores outside of HLA region using idx from foverlaps
    ans <- all.scoresinfo[!idx[!is.na(yid), xid]]
}

#' In case of immune cells, we do not perform any further annotation and just
#' sum up all scores to generate a single genome-wide score also including cis
#' effects.
##------------------------------------------------------------------------------
if (analysis == "immunecells")) {
    cat("Immune cells cannot be currently annotated. Skipping.\n")
    ans[, trait_name := gsub(":", "_", trait_name)]
    ans[, trait_name := gsub(" ", "_", trait_name)]

    ## append X_ to scorenames
    ans[, scoreid := paste0("X_", scoreid)]
    scorenames <- colnames(all.scores)
    x.scorenames <- paste0("X_", scorenames)
    colnames(all.scores) <- x.scorenames

    ## keep all cells which have only one score
    n.scores <- ans[minpvalue<1E-6, .N, by=c("trait_name", "gwasid")]
    cells.with.single.score <- n.scores[N==1]
    single.scores <- ans[cells.with.single.score,
        on=c("gwasid", "trait_name")][minpvalue<1E-6]
    scores.matrix <- all.scores[, colnames(all.scores) %in% single.scores$scoreid]
    colnames(scores.matrix) <- single.scores$trait_name

    ## sum all scores for each cell trait
    cell.traits <- n.scores[N>1, trait_name]
    genome.wide.scores <- NULL
    for (this.trait in cell.traits) {
        this.scoreids <- ans[trait_name %in% this.trait, scoreid]
        this.scores <- matrix(rowSums(
            all.scores[, colnames(all.scores) %in% this.scoreids]))
        colnames(this.scores) <- this.trait
        genome.wide.scores <- cbind(genome.wide.scores, this.scores)
    }
    genome.wide.scores <- cbind(genome.wide.scores, scores.matrix)

    ## save processed scores
    save(genome.wide.scores,
         file=file.path(output.dir, "trans.genotypicscore.1e-6.Rdata.gz"))
    save(genome.wide.scores,
         file=file.path(output.dir, "trans.scoresinfo.1e-6.Rdata.gz"))
}

#' eQTL and pQTL scores have to be annotated and classified into cis, cis-x and
#' trans. We have a function in GENOSCORES package to annotate eQTLs which can
#' be called directly. Annotation of pQTLs can currently be performed only on
#' GENOSCORES server, therefore, this functionality is disabled. Contact us if
#' you need help with annotating the pQTLs.
##------------------------------------------------------------------------------
if (analysis == "eQTL") ans <- annotate.qtls(ans)

if (analysis == "pQTL")
    stop("pQTLs cannot be annotated via API. Contact us for help.\n")

## DO NOT RUN. This block can be used only on GENOSCORES server directly. You
## can inspect `annotate.pqtls.R` and write a similar script to manually
## annotate pQTL scores if you have a uniprot table available locally.
if (FALSE) {
    ## TODO: create API function to annotate pQTLs (see !217)
    source("annotate.pqtls.R")
    genoscores.connect(3308)
    ans <- annotate.scoresinfo(ans)
    genoscores.disconnect()
}

#' Aggregate annotated scores into genome-wide trans scores.
##------------------------------------------------------------------------------

## create list of candidates: get all genes with trans signals
trans.genes <- ans[qtl_type == "trans", unique(gene_symbol)]

## remove scores for genes not in a candidate list
ans <- ans[gene_symbol %in% trans.genes]

## append X_ to scorenames
ans[, scoreid := paste0("X_", scoreid)]
scorenames <- colnames(all.scores)
x.scorenames <- paste0("X_", scorenames)
colnames(all.scores) <- x.scorenames

## filter scores matrix
fs <- all.scores[, colnames(all.scores) %in% ans[, scoreid]]

## compute variance of each locus-specific score
score.var <- apply(fs, 2, var)
score.var <- setDT(stack(score.var))
setnames(score.var, old=c("values", "ind"), new=c("variance", "scoreid"))

## merge with scoresinfo
ans <- ans[score.var, on="scoreid", nomatch=NULL]

## order scoresinfo to respect scores matrix
ans <- ans[match(colnames(fs), ans[, scoreid]), ]

## check that scoresinfo is identical to scores matrix
stopifnot(identical(colnames(fs), ans[, scoreid]))

## group by gene
ans <- ans[, .SD, by=gene_symbol]

## sum trans scores and keep cis score corresponding to each gene
genome.wide.scores <- NULL
cis.scoreids <- ans[qtl_type=="cis" | qtl_type=="cis-x", scoreid]
genome.wide.scores <- fs[, colnames(fs) %in% cis.scoreids]
## we proceed by gwasid computing genome-wide trans scores per GWAS
## TODO: see issue #27
n.trans <- ans[qtl_type=="trans" & minpvalue<1E-6, .N, by="gwasid"]
gwas.with.single.score <- n.trans[N==1, gwasid]
single.trans.scoreids <- ans[qtl_type=="trans" & minpvalue<1E-6 &
                             gwasid %in% gwas.with.single.score, scoreid]
trans.matrix <- fs[, colnames(fs) %in% single.trans.scoreids]

## create a trans score id as X_<gwasid>_trans
ans[scoreid %in% single.trans.scoreids, scoreid := sprintf("X_%s_trans", gwasid)]
scorenames <- ans[grep("trans", scoreid), scoreid]
colnames(trans.matrix) <- scorenames
genome.wide.scores <- cbind(genome.wide.scores, trans.matrix)

## sanity checks
stopifnot(all(colnames(trans.matrix) %in% colnames(genome.wide.scores)))
stopifnot(all(colnames(trans.matrix) %in% ans[, scoreid]))

## sum up multiple trans scores per gwas
gwas.with.many.scores <- n.trans[N>1]
for (this.gwasid in gwas.with.many.scores[, gwasid]) {
    this.scoreids <- ans[qtl_type=="trans" & minpvalue<1E-6 &
                         gwasid==this.gwasid, scoreid]
    this.combined.score <- matrix(rowSums(fs[, colnames(fs) %in% this.scoreids]))
    ans[scoreid %in% this.scoreids, scoreid := sprintf("X_%s_trans", this.gwasid)]
    colnames(this.combined.score) <- sprintf("X_%s_trans", this.gwasid)
    genome.wide.scores <- cbind(genome.wide.scores, this.combined.score)
}

ans <- ans[order(gene_symbol)]

#' Check that all scores were aggregated correctly.
##------------------------------------------------------------------------------
names <- colnames(genome.wide.scores)
cis <- names[grep("X_.*_.*_.*", names)]
trans <- names[grep("X_.*_trans", names)]
stopifnot(identical(cis.scoreids[order(cis.scoreids)], cis[order(cis)]))
single.trans <- colnames(trans.matrix)
multiple.trans <- sprintf("X_%s_trans", gwas.with.many.scores$gwasid)
all.trans <- c(single.trans, multiple.trans)
stopifnot(identical(all.trans[order(all.trans)], trans[order(trans)]))

#' Save processed scores.
##------------------------------------------------------------------------------
trans.genome.wide.scoresinfo <- ans
save(genome.wide.scores,
     file=file.path(output.dir, "trans.genotypicscore.1e-6.Rdata.gz"))
save(trans.genome.wide.scoresinfo,
     file=file.path(output.dir, "trans.scoresinfo.1e-6.Rdata.gz"))
