#' This is an example script which loads either eQTL or pQTL scores depending
#' on the value of `analysis` variable. In addition, this script requires that
#' the path to the scores directory is specified in the global environment.
## =============================================================================
if (analysis == "eQTL") {
    ## load processed scores and metadata
    eqtl.datadir <- score.dir
    load(file.path(eqtl.datadir, "trans.scoresinfo.1e-6.Rdata.gz"))
    load(file.path(eqtl.datadir, "trans.genotypicscore.1e-6.Rdata.gz"))

    ## clean up metadata
    trans.genome.wide.scoresinfo[, chrom_min :=
        factor(chrom_min, levels=c(1:22, "X", "Y"))]

    trans.genome.wide.scoresinfo[, x.startpos :=
        position.absolute(chrom_min, startpos) - DMhits.flanksize.kb * 1000]
    trans.genome.wide.scoresinfo[, x.endpos :=
        position.absolute(chrom_min, endpos) + DMhits.flanksize.kb * 1000]
    setkey(trans.genome.wide.scoresinfo, x.startpos, x.endpos)
    setkey(t1d.hits, x.clumpStart, x.clumpEnd)

    eqtl.overlaps <- foverlaps(trans.genome.wide.scoresinfo[, .(scoreid, x.startpos, x.endpos)],
                               t1d.hits[, .(nearest, x.clumpStart, x.clumpEnd)], nomatch=NULL)
    eqtl.overlaps <- eqtl.overlaps[, .(near.t1dgenes=paste(nearest, collapse=",")), by=scoreid]

    setkey(eqtl.overlaps, scoreid)
    setkey(trans.genome.wide.scoresinfo, scoreid)
    trans.genome.wide.scoresinfo <- eqtl.overlaps[trans.genome.wide.scoresinfo]

    setorder(trans.genome.wide.scoresinfo, chrom_min, startpos, endpos)
    transqtl.loci <- trans.genome.wide.scoresinfo[qtl_type=="trans",
        .(scoreid,
          gene_symbol,
          chrom_min,
          startpos,
          endpos,
          nearby_genes,
          near.t1dgenes)]
    transqtl.loci[, x.startpos := position.absolute(chrom_min, startpos)]
    transqtl.loci[, x.endpos := position.absolute(chrom_min, endpos)]
    setorder(transqtl.loci, chrom_min, x.startpos)

    transqtl.loci[grep("^STH|NSF$", nearby_genes)]
    ## one score clump has STH included as nearby.t1dgene though it is >600 kb downstream of transcription site
    ## this is because SNPs in the T1D GWAS hit extend from 45.32 to 46.50 Mb

    ## define region by minimum gap of between last endpos and startpos
    transqtl.loci[, newregion :=
        c(1, as.integer((x.startpos[-1] - x.endpos[-.N] > 1E5)))]
    transqtl.loci[, region := cumsum(newregion)]
    setkey(transqtl.loci, scoreid)
    setkey(trans.genome.wide.scoresinfo, scoreid)
    trans.genome.wide.scoresinfo <-
        transqtl.loci[, .(scoreid, region)][trans.genome.wide.scoresinfo]
    setorder(transqtl.loci, chrom_min, startpos)

    ## collapse to one row per region as transqtl.regions
    transqtl.regions <- transqtl.loci[,
        .(chrom=chrom_min[1],
          startpos=min(startpos),
          endpos=max(endpos),
          numscores=.N,
          targetgenes=paste(unique(gene_symbol), collapse=","),
          nearbygenes=paste(unique(unlist(strsplit(nearby_genes, split=","))), collapse=","),
          near.t1dgenes = paste(unique(unlist(strsplit(near.t1dgenes, split=","))), collapse = ",")),
        by=region]

    setkey(transqtl.regions, region)
    transqtl.regions[nchar(near.t1dgenes) > 0] # regions containing known T1D-associated SNPs

    ## generate metadata table with rows matching colnames of scores matrix:
    ## i.e. genomewide trans scores
    ## for genomewide trans scores the scores matrix colname is the gene symbol
    ## for cis scores the scores matrix colname is the locus-specific scoreid
    tscores.info <- trans.genome.wide.scoresinfo[qtl_type=="trans",
        .(numscores=.N,
          matrix.colname=paste0("X_",gwasid,"_trans"),
          qtl_type="trans",
          gene_symbol=gene_symbol[1],
          gene_chrom=gene_chrom[1],
          gene_startpos=gene_startpos[1],
          gene_endpos=gene_endpos[1],
          regions=paste(sort(unique(region)),
                        collapse=","),
          locus.diversity=diversity(variance, 1)),
        by=.(traitid, gwasid)]

    cscores.info <- trans.genome.wide.scoresinfo[substr(qtl_type, 1, 3)=="cis",
        .(scoreid,
          matrix.colname=scoreid,
          qtl_type,
          gene_symbol,
          gene_chrom,
          gene_startpos, gene_endpos)]

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
}

if (analysis == "pQTL") {
    ## load processed scores and metadata
    pqtl.datadir <- score.dir
    load(file.path(pqtl.datadir, "trans.scoresinfo.1e-6.Rdata.gz"))
    load(file.path(pqtl.datadir, "trans.genotypicscore.1e-6.Rdata.gz"))

    ## rename and clean up metadata
    ## cis- columns have name X_nnn_n_n
    ## trans- columns  have gene name
    trans.genome.wide.scoresinfo[, gene_symbol := gsub("character\\(0\\),?", "", gene_symbol)]
    trans.genome.wide.scoresinfo[, gene_symbol := gsub("c\\(", "", gene_symbol)]
    trans.genome.wide.scoresinfo[, gene_symbol := gsub("\\)", "", gene_symbol)]
    trans.genome.wide.scoresinfo[, gene_symbol := gsub('\\"', "", gene_symbol)]
    trans.genome.wide.scoresinfo[, gene_symbol := gsub(" TBCE", "", gene_symbol)]
    trans.genome.wide.scoresinfo[, gene_symbol := gsub(",$", "", gene_symbol)]
    trans.genome.wide.scoresinfo[is.na(gene_chrom), gene_symbol]
    trans.genome.wide.scoresinfo <- trans.genome.wide.scoresinfo[nchar(gene_symbol) > 0]

    ## FIXME: some trans scores do not map to a single gene
    trans.genome.wide.scoresinfo[, gene_symbol := gsub(",.*", "", gene_symbol)]
    pqtl.genes.missinginfo <- trans.genome.wide.scoresinfo[is.na(gene_chrom), gene_symbol]
    missinginfo <- genes.info[match(pqtl.genes.missinginfo, gene_symbol),
                                       .(chromosome, gene_startpos, gene_endpos)]

    trans.genome.wide.scoresinfo[is.na(gene_chrom),
        `:=`(gene_chrom=missinginfo$chromosome,
             gene_startpos=missinginfo$gene_startpos,
             gene_endpos=missinginfo$gene_endpos)]

    options(warn=-1)
    trans.genome.wide.scoresinfo[, gene_chrom :=
        factor(as.character(gene_chrom), levels=c(1:22, "X", "Y"))]
    trans.genome.wide.scoresinfo[, gene_startpos := as.integer(gene_startpos)]
    trans.genome.wide.scoresinfo[, gene_endpos := as.integer(gene_endpos)]
    trans.genome.wide.scoresinfo[, chrom_min :=
        factor(chrom_min, levels=c(1:22, "X", "Y"))]
    options(warn=2)

    summary(trans.genome.wide.scoresinfo)

    ## trans scores identified in matrix only by gene_symbol, so may have been
    ## summed over multiple gwasids
	pqtl.tscores.info <- unique(trans.genome.wide.scoresinfo[qtl_type=="trans",
        .(numscores=.N,
          scoreids=paste(scoreid, collapse=","),
          gwasids=paste(gwasid, collapse=","),
          numgwas=length(unique(gwasid)),
          matrix.colname=paste0("X_",gwasid,"_trans"),
          qtl_type="trans",
          gene_symbol=gene_symbol[1],
          gene_chrom=gene_chrom[1],
          gene_startpos=gene_startpos[1],
          gene_endpos=gene_endpos[1],
          regions=paste(sort(unique(region)),
                        collapse=","),
          locus.diversity=diversity(variance, 1)),
        by=gene_symbol])

    ## cis scores are uniquely identified by scoreid
    pqtl.cscores.info <- trans.genome.wide.scoresinfo[substr(qtl_type, 1, 3)=="cis",
        .(matrix.colname=scoreid,
          scoreids=paste(scoreid, collapse=","),
          gwasids=gwasid,
          qtl_type,
          gene_symbol,
          gene_chrom,
          gene_startpos, gene_endpos),
        by=.(traitid, gwasid)]

    pqtl.allscores.info <- rbind(pqtl.tscores.info, pqtl.cscores.info, fill=TRUE)
}
