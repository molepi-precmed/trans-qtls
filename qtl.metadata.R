#' This is an example script which loads either eQTL or pQTL scores depending
#' on the value of `analysis` variable. In addition, this script requires that
#' the path to the scores directory is specified in the global environment.
## =============================================================================
if (analysis == "eQTL") {
    ## load processed scores and metadata
    eqtl.datadir <- score.dir
    load(paste0(eqtl.datadir, "trans.scoresinfo.1e-6.Rdata.gz"))
    load(paste0(eqtl.datadir, "trans.genotypicscore.1e-6.Rdata.gz"))

    ## rename and clean up metadata
    eqtl.info <- trans.genome.wide.scoresinfo
    eqtl.info[, chrom_min := factor(chrom_min, levels=c(1:22, "X", "Y"))]

    eqtl.info[, x.startpos := position.absolute(chrom_min, startpos) - DMhits.flanksize.kb * 1000]
    eqtl.info[, x.endpos := position.absolute(chrom_min, endpos) + DMhits.flanksize.kb * 1000]
    setkey(eqtl.info, x.startpos, x.endpos)
    setkey(t1d.hits, x.clumpStart, x.clumpEnd)

    eqtl.overlaps <- foverlaps(eqtl.info[, .(scoreid, x.startpos, x.endpos)],
                               t1d.hits[, .(nearest, x.clumpStart, x.clumpEnd)], nomatch=NULL)
    eqtl.overlaps <- eqtl.overlaps[, .(near.t1dgenes=paste(nearest, collapse=",")), by=scoreid]

    setkey(eqtl.overlaps, scoreid)
    setkey(eqtl.info, scoreid)
    eqtl.info <- eqtl.overlaps[eqtl.info]

    ## transqtl.loci has 15887 loci contributing to trans- scores
    setorder(eqtl.info, chrom_min, startpos, endpos)
    transqtl.loci <- eqtl.info[qtl_type=="trans", .(scoreid, gene_symbol, chrom_min, startpos, endpos, nearby_genes, near.t1dgenes)]
    transqtl.loci[, x.startpos := position.absolute(chrom_min, startpos)]
    transqtl.loci[, x.endpos := position.absolute(chrom_min, endpos)]
    setorder(transqtl.loci, chrom_min, x.startpos)

    transqtl.loci[grep("^STH|NSF$", nearby_genes)]
    ## one score clump has STH included as nearby.t1dgene though it is >600 kb downstream of transcription site
    ## this is because SNPs in the T1D GWAS hit extend from 45.32 to 46.50 Mb

    ## X_632743_1_17       ERMAP        17 46678191 46681485  NSF   STH,NSF 2537458753 2537462047

    ## define region by minimum gap of between last endpos and startpos
    transqtl.loci[, newregion := c(1, as.integer((x.startpos[-1] - x.endpos[-.N] > 1E5)))]
    transqtl.loci[, region := cumsum(newregion)]
    setkey(transqtl.loci, scoreid)
    setkey(eqtl.info, scoreid)
    eqtl.info <-
        transqtl.loci[, .(scoreid, region)][eqtl.info]
    setorder(transqtl.loci, chrom_min, startpos)

    ## collapse to one row per region as transqtl.regions (1198 regions)
    transqtl.regions <- transqtl.loci[, .(chrom=chrom_min[1],
                                          startpos=min(startpos),
                                          endpos=max(endpos),
                                          numscores=.N,
                                          targetgenes=paste(unique(gene_symbol), collapse=","),
                                          nearbygenes=paste(unique(unlist(strsplit(nearby_genes, split=","))), collapse=","),
                                          near.t1dgenes = paste(unique(unlist(strsplit(near.t1dgenes, split=","))), collapse = ",")
                                          ),
                                      by=region]

    setkey(transqtl.regions, region)
    transqtl.regions[nchar(near.t1dgenes) > 0] # regions containing known T1D-associated SNPs

    ## generate metadata table with rows matching colnames of scores matrix: i.e. genomewide trans scores
    ## for genomewide trans scores the scores matrix colname is the gene symbol
    ## for cis scores the scores matrix colname is the locus-specific scoreid
    tscores.info <- eqtl.info[qtl_type=="trans",
                                                 .(numscores=.N,
                                                   matrix.colname=gene_symbol[1],
                                                   qtl_type="trans",
                                                   gene_symbol=gene_symbol[1],
                                                   gene_chrom=gene_chrom[1],
                                                   gene_startpos=gene_startpos[1],
                                                   gene_endpos=gene_endpos[1],
                                                   regions=paste(sort(unique(region)),
                                                                 collapse=","),
                                                   locus.diversity=diversity(variance, 1)),
                                                 by=.(traitid, gwasid)]
    assign("tscores.info", tscores.info, envir=.GlobalEnv)
    cscores.info <- eqtl.info[substr(qtl_type, 1, 3)=="cis",
                                                 .(scoreid,
                                                   matrix.colname=scoreid,
                                                   qtl_type,
                                                   gene_symbol,
                                                   gene_chrom,
                                                   gene_startpos, gene_endpos)]
    assign("cscores.info", cscores.info, envir=.GlobalEnv)

    ## combine the trans and cis scores
    allscores.info <- rbind(tscores.info, cscores.info, fill=TRUE)
    ## column matrix.colname in allscores.info matches column names in scores table
    ## recode cis-x as cis but keep an indicator for cis-x
    allscores.info[, cisx := qtl_type=="cis-x"]
    allscores.info[qtl_type=="cis-x", qtl_type := "cis"]
    setnames(allscores.info, "gene_chrom", "chromosome")

    table(allscores.info$matrix.colname %in% colnames(genome.wide.scores))
    ## remove rows not matched in scores table
    allscores.info <- allscores.info[matrix.colname %in% colnames(genome.wide.scores)]

    scores.sd <- apply(genome.wide.scores, 2, sd)
    scores.sd <- data.table(matrix.colname=names(scores.sd), sdscore=as.numeric(scores.sd))
    setkey(scores.sd, matrix.colname)
    setkey(allscores.info, matrix.colname)
    allscores.info <- scores.sd[allscores.info]
    assign("allscores.info", allscores.info, envir=.GlobalEnv)
}

if (analysis == "pQTL") {
    ## load processed scores and metadata
    pqtl.dir <- score.dir
    load(paste0(pqtl.datadir, "trans.scoresinfo.1e-6.Rdata.gz"))
    load(paste0(pqtl.datadir, "trans.genotypicscore.1e-6.Rdata.gz"))

    pqtl.scores <- genome.wide.scores # 3849 scores
    assign("pqtl.scores", pqtl.scores, envir=.GlobalEnv)

    ## rename and clean up metadata
    ## cis- columns have name X_nnn_n_n
    ## trans- columns  have gene name
    pqtl.info <- trans.genome.wide.scoresinfo # scoreid
    pqtl.info[, gene_symbol := gsub("character\\(0\\),?", "", gene_symbol)]
    pqtl.info[, gene_symbol := gsub("c\\(", "", gene_symbol)]
    pqtl.info[, gene_symbol := gsub("\\)", "", gene_symbol)]
    pqtl.info[, gene_symbol := gsub('\\"', "", gene_symbol)]
    pqtl.info[, gene_symbol := gsub(" TBCE", "", gene_symbol)]
    pqtl.info[, gene_symbol := gsub(",$", "", gene_symbol)]
    pqtl.info[is.na(gene_chrom), gene_symbol]
    pqtl.info <- pqtl.info[nchar(gene_symbol) > 0]

    ##FIXME: some trans scores do not map to a single gene
    pqtl.info[, gene_symbol := gsub(",.*", "", gene_symbol)]
    pqtl.genes.missinginfo <- pqtl.info[is.na(gene_chrom), gene_symbol]
    missinginfo <- genes.info[match(pqtl.genes.missinginfo, gene_symbol),
                                       .(chromosome, gene_startpos, gene_endpos)]

    pqtl.info[is.na(gene_chrom), `:=`(gene_chrom=missinginfo$chromosome,
                                      gene_startpos=missinginfo$gene_startpos,
                                      gene_endpos=missinginfo$gene_endpos)]

    options(warn=-1)
    pqtl.info[, gene_chrom := factor(as.character(gene_chrom), levels=c(1:22, "X", "Y"))]
    pqtl.info[, gene_startpos := as.integer(gene_startpos)]
    pqtl.info[, gene_endpos := as.integer(gene_endpos)]
    pqtl.info[, chrom_min := factor(chrom_min, levels=c(1:22, "X", "Y"))]
    options(warn=2)

    # 3049 unique traitids
    summary(pqtl.info) # 56072 trans- scores, 793 cis, 302 cis-x

    ## trans scores identified in matrix only by gene_symbol, so may have been summed over multiple gwasids
    # 2983 trans scores
    pqtl.tscores.info <- unique(pqtl.info[qtl_type=="trans",
                                          .(numscores=.N,
                                            scoreids=paste(scoreid, collapse=","),
                                            gwasids=paste(gwasid, collapse=","),
                                            numgwas=length(unique(gwasid)),
                                            matrix.colname=gene_symbol[1],
                                            qtl_type="trans",
                                            gene_symbol=gene_symbol[1],
                                            gene_chrom=gene_chrom[1],
                                            gene_startpos=gene_startpos[1],
                                            gene_endpos=gene_endpos[1],
                                            regions=paste(sort(unique(region)),
                                                          collapse=","),
                                            locus.diversity=diversity(variance, 1)),
                                          by=gene_symbol])
    assign("pqtl.tscores.info", pqtl.tscores.info, envir=.GlobalEnv)

    ## cis scores are uniquely identified by scoreid
    pqtl.cscores.info <- pqtl.info[substr(qtl_type, 1, 3)=="cis",
                                   .(matrix.colname=scoreid,
                                     scoreids=paste(scoreid, collapse=","),
                                     gwasids=gwasid,
                                     qtl_type,
                                     gene_symbol,
                                     gene_chrom,
                                     gene_startpos, gene_endpos),
                                   by=.(traitid, gwasid)]
    assign("pqtl.cscores.info", pqtl.cscores.info, envir=.GlobalEnv)

    pqtl.allscores.info <- rbind(pqtl.tscores.info, pqtl.cscores.info, fill=TRUE)
    assign("pqtl.allscores.info", pqtl.allscores.info, envir=.GlobalEnv)
}
