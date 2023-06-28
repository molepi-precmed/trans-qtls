#' This is an example script which loads either eQTL or pQTL scores depending
#' on the value of `analysis` variable. In addition, this script requires that
#' the path to the scores directory is specified in the global environment.
## =============================================================================
if (analysis == "eQTL") {
  ## load processed scores and metadata
  eqtl.datadir <- score.dir
  load(file.path(eqtl.datadir, "trans.scoresinfo.1e-6.Rdata.gz"))
  load(file.path(eqtl.datadir, "trans.genotypicscore.1e-6.Rdata.gz"))
  load(file.path(eqtl.datadir, "all.scoresinfo.annotated.Rdata.gz"))
  
  ## clean up metadata
  trans.genome.wide.scoresinfo[, chrom_min :=
                                 factor(chrom_min, levels=c(1:22, "X", "Y"))]
  
  trans.genome.wide.scoresinfo[, x.startpos :=
                                 position.absolute(chrom_min, startpos) - DMhits.flanksize.kb * 1000]
  trans.genome.wide.scoresinfo[, x.endpos :=
                                 position.absolute(chrom_min, endpos) + DMhits.flanksize.kb * 1000]
  setkey(trans.genome.wide.scoresinfo, x.startpos, x.endpos)
  setkey(t1d.hits, x.clumpStart, x.clumpEnd)
  
  eqtl.overlaps <- foverlaps(
    trans.genome.wide.scoresinfo[, .(scoreid, x.startpos, x.endpos)],
    t1d.hits[, .(nearest, x.clumpStart, x.clumpEnd)], nomatch=NULL)
  eqtl.overlaps <- eqtl.overlaps[,
                                 .(near.t1dgenes=paste(nearest, collapse=",")), by=scoreid]
  
  setkey(eqtl.overlaps, scoreid)
  setkey(trans.genome.wide.scoresinfo, scoreid)
  trans.genome.wide.scoresinfo <- eqtl.overlaps[trans.genome.wide.scoresinfo]
  
  ## join with all annotated score info to bring in contributing scoreids
  ans[, x.startpos :=
        position.absolute(chrom_min, startpos) - DMhits.flanksize.kb * 1000]
  ans[, x.endpos :=
        position.absolute(chrom_min, endpos) + DMhits.flanksize.kb * 1000]
  setkey(ans, x.startpos, x.endpos)
  
  trans.genome.wide.scoresinfo[ans[, .(gene_symbol, x.startpos, x.endpos,
                                       scoreid)],
                               on=c("gene_symbol", "x.startpos", "x.endpos"),
                               contributing.score := i.scoreid]
  eqtl.trans.genome.wide.scoresinfo <- trans.genome.wide.scoresinfo
  
  setorder(eqtl.trans.genome.wide.scoresinfo, chrom_min, startpos, endpos)
  
  transqtl.loci <- eqtl.trans.genome.wide.scoresinfo[
    qtl_type=="trans",
    .(scoreid,
      gene_symbol,
      chrom_min,
      startpos,
      endpos,
      x.startpos,
      x.endpos,
      nearby_genes,
      near.t1dgenes)]
  setorder(transqtl.loci, chrom_min, x.startpos)
  
  ## define region by minimum gap of 1E5 between last endpos and startpos
  transqtl.loci[, newregion :=
                  c(1, as.integer((x.startpos[-1] - x.endpos[-.N] > 1E5)))]
  transqtl.loci[, region := cumsum(newregion)]
  
  ## merge region field from transqtl.loci with eqtl.info, by scoreid
  eqtl.trans.genome.wide.scoresinfo[transqtl.loci,
                                    on=c("scoreid", "x.startpos", "x.endpos"),
                                    region := i.region]
  setorder(transqtl.loci, chrom_min, x.startpos)
  
  ## collapse to one row per region as transqtl.regions
  transqtl.regions <- transqtl.loci[,
                                    .(chrom=chrom_min[1],
                                      startpos=min(startpos),
                                      endpos=max(endpos),
                                      numscores=.N,
                                      targetgenes=paste(unique(gene_symbol), collapse=","),
                                      nearbygenes=paste(unique(unlist(strsplit(nearby_genes, split=","))), collapse=","),
                                      near.t1dgenes=paste(unique(unlist(strsplit(near.t1dgenes, split=","))), collapse = ",")),
                                    by=region] ## 602
  
  setkey(transqtl.regions, region)
  transqtl.regions[nchar(near.t1dgenes) > 0] # regions containing known SLE-associated SNPs
  
  ## generate metadata table with rows matching colnames of scores matrix:
  ## i.e. genomewide trans scores
  ## for genomewide trans scores the scores matrix colname is the gene symbol
  ## for cis scores the scores matrix colname is the locus-specific scoreid
  tscores.info <- eqtl.trans.genome.wide.scoresinfo[
    qtl_type=="trans",
    .(numscores=.N,
      matrix.colname=scoreid,
      qtl_type="trans",
      gene_symbol=gene_symbol[1],
      gene_chrom=gene_chrom[1],
      gene_startpos=gene_startpos[1],
      gene_endpos=gene_endpos[1],
      regions=paste(sort(unique(region)),
                    collapse=","),
      locus.diversity=diversity(variance, 1)),
    by=.(traitid, gwasid)]
  
  cscores.info <- eqtl.trans.genome.wide.scoresinfo[
    substr(qtl_type, 1, 3)=="cis",
    .(scoreid,
      matrix.colname=scoreid,
      qtl_type,
      gene_symbol,
      gene_chrom,
      gene_startpos, gene_endpos)]
  
  ## combine the trans and cis scores
  allscores.info <- rbind(tscores.info, cscores.info, fill=TRUE)
  ## column matrix.colname in allscores.info matches column names in scores
  ## table. Recode cis-x as cis but keep an indicator for cis-x
  allscores.info[, cisx := qtl_type=="cis-x"]
  allscores.info[qtl_type=="cis-x", qtl_type := "cis"]
  
  table(allscores.info$matrix.colname %in% colnames(genome.wide.scores))
  ## remove rows not matched in scores table (all matched)
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
  load(file.path(pqtl.datadir, "all.scoresinfo.annotated.Rdata.gz"))
  
  ## rename and clean up metadata
  ## cis- columns have name X_nnn_n_n
  ## trans- columns  have gene name
  trans.genome.wide.scoresinfo[, gene_symbol := gsub("character\\(0\\),?", "", gene_symbol)]
  trans.genome.wide.scoresinfo[, gene_symbol := gsub("c\\(", "", gene_symbol)]
  trans.genome.wide.scoresinfo[, gene_symbol := gsub("\\)", "", gene_symbol)]
  trans.genome.wide.scoresinfo[, gene_symbol := gsub('\\"', "", gene_symbol)]
  trans.genome.wide.scoresinfo[, gene_symbol := gsub(" TBCE", "", gene_symbol)]
  trans.genome.wide.scoresinfo[, gene_symbol := gsub(",$", "", gene_symbol)]
  trans.genome.wide.scoresinfo[gene_chrom=="c(\"1\", \"1\")", gene_chrom := 1]
  trans.genome.wide.scoresinfo[gene_chrom=="c(\"19\", \"19\")", gene_chrom := 19]
  trans.genome.wide.scoresinfo[, gene_startpos := gsub("c\\((\\d*),.*\\)", "\\1", gene_startpos)]
  trans.genome.wide.scoresinfo[, gene_endpos := gsub("c\\((\\d*),.*\\)", "\\1", gene_endpos)]
  trans.genome.wide.scoresinfo[, gene_symbol := gsub("(.*),.*", "\\1", gene_symbol)]
  trans.genome.wide.scoresinfo[is.na(gene_chrom), gene_symbol]
  
  trans.genome.wide.scoresinfo[, gene_chrom :=
                                 factor(gene_chrom, levels=c(1:22, "X", "Y"))]
  trans.genome.wide.scoresinfo[, gene_startpos := as.integer(gene_startpos)]
  trans.genome.wide.scoresinfo[, gene_endpos := as.integer(gene_endpos)]
  
  trans.genome.wide.scoresinfo[, x.startpos := position.absolute(chrom_min, startpos) - DMhits.flanksize.kb * 1000]
  trans.genome.wide.scoresinfo[, x.endpos := position.absolute(chrom_min, endpos) + DMhits.flanksize.kb * 1000]
  
  ## FIXME: some trans scores do not map to a single gene
  trans.genome.wide.scoresinfo[, gene_symbol := gsub(",.*", "", gene_symbol)]
  
  trans.genome.wide.scoresinfo[, gene_chrom :=
                                 factor(as.character(gene_chrom), levels=c(1:22, "X", "Y"))]
  trans.genome.wide.scoresinfo[, chrom_min :=
                                 factor(chrom_min, levels=c(1:22, "X", "Y"))]
  
  summary(trans.genome.wide.scoresinfo)
  
  ## crate table of overlaps with known T1D genes
  setkey(trans.genome.wide.scoresinfo, chrom_min, x.startpos, x.endpos)
  setkey(t1d.hits, chromosome, x.clumpStart, x.clumpEnd)
  pqtl.overlaps <- foverlaps(
    trans.genome.wide.scoresinfo[, .(scoreid, chrom_min, x.startpos, x.endpos)],
    t1d.hits[, .(nearest, chromosome, x.clumpStart, x.clumpEnd)],
    nomatch=NULL)
  
  pqtl.overlaps <- pqtl.overlaps[,
                                 .(near.t1dgenes=paste(unique(nearest), collapse=",")),
                                 by=c("chrom_min", "x.startpos", "x.endpos")]
  
  ## merge back with metadata
  pqt.info <- pqtl.overlaps[trans.genome.wide.scoresinfo,
                            on=c("chrom_min", "x.startpos", "x.endpos")]
  pqt.info[, region := NULL]
  
  ## join with all annotated score info to bring in contributing scoreids
  ans[, x.startpos := position.absolute(chrom_min, startpos)]
  ans[, x.endpos := position.absolute(chrom_min, endpos)]
  pqt.info[ans[, .(gene_symbol, x.startpos, x.endpos, scoreid)],
           on=c("gene_symbol", "x.startpos", "x.endpos"),
           contributing.score := i.scoreid]
  
  ## create transpqtl.loci
  setorder(pqt.info, chrom_min, startpos, endpos)
  transpqtl.loci <- pqt.info[qtl_type=="trans",
                             .(scoreid, gene_symbol, chrom_min, startpos,
                               endpos, x.startpos, x.endpos, nearby_genes,
                               near.t1dgenes, studyid)]
  setorder(transpqtl.loci, chrom_min, x.startpos)
  
  ## define regions by minimum gap of 1E5 between last endpos and startpos
  transpqtl.loci[, newregion :=
                   c(1, as.integer((x.startpos[-1] - x.endpos[-.N] > 1E5)))]
  transpqtl.loci[, region := cumsum(newregion)]
  
  ## merge region field from transqtl.loci with eqtl.info, by scoreid
  pqt.info[transpqtl.loci, on=c("scoreid", "x.startpos", "x.endpos"),
           region := i.region]
  setorder(transpqtl.loci, chrom_min, x.startpos)
  
  ## collapse to one row per region as transpqtl.regions
  transpqtl.regions <- transpqtl.loci[,
                                      .(chrom=chrom_min[1],
                                        startpos=min(startpos),
                                        endpos=max(endpos),
                                        numscores=.N,
                                        targetgenes=paste(unique(gene_symbol), collapse=","),
                                        nearbygenes=paste(unique(unlist(strsplit(nearby_genes, split=","))), collapse=","),
                                        near.t1dgenes=paste(unique(unlist(strsplit(near.t1dgenes, split=","))), collapse = ",")),
                                      by=region] ## 1086
  
  ## generate metadata table with rows matching colnames of scores matrix:
  ## i.e. genomewide trans scores
  ## for genomewide trans scores the scores matrix colname is the gene symbol
  ## for cis scores the scores matrix colname is the locus-specific scoreid
  pqtl.tscores.info <- unique(pqt.info[
    qtl_type=="trans",
    .(numscores=.N,
      scoreids=paste(scoreid, collapse=","),
      gwasids=paste(gwasid, collapse=","),
      numgwas=length(unique(gwasid)),
      matrix.colname=paste0("X_", gwasid, "_trans"),
      qtl_type,
      trait_name=trait_name[1],
      gene_symbol=gene_symbol[1],
      gene_chrom=gene_chrom[1],
      gene_startpos=gene_startpos[1],
      gene_endpos=gene_endpos[1],
      regions=paste(sort(unique(region)),
                    collapse=","),
      locus.diversity=diversity(variance, 1)),
    by=gwasid])
  
  ## cis scores are uniquely identified by scoreid
  pqtl.cscores.info <- pqt.info[
    substr(qtl_type, 1, 3)=="cis",
    .(matrix.colname=scoreid,
      scoreids=paste(scoreid, collapse=","),
      studyid,
      gwasid,
      qtl_type,
      trait_name,
      gene_symbol,
      gene_chrom,
      gene_startpos,
      gene_endpos),
    by=gwasid]
  
  pqtl.allscores.info <- rbind(pqtl.tscores.info, pqtl.cscores.info, fill=TRUE)
  
  scores.sd <- apply(genome.wide.scores, 2, sd)
  scores.sd <- data.table(matrix.colname=names(scores.sd), sdscore=as.numeric(scores.sd))
  setkey(scores.sd, matrix.colname)
  setkey(allscores.info, matrix.colname)
  pqtl.allscores.info <- scores.sd[pqtl.allscores.info]
}
