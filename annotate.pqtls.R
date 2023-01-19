## This is a blueprint of a function to annotate protein scores.
## It uses proteins2genes() function to get all protein annotations from
## GENOSCORES database and then intelligently intersects these annotations with
## the scoresinfo object provided by the user.
##
## Author: Athina Spiliopoulou, Andrii Iakovliev
## Date: 7 July 2022
## =============================================================================

#' Concatenate values from multiple rows of uniprotinfo, if a protein maps to
#' multiple uniprot ids
#'
#' @param ids.list vector with mapped uniprot-ids for a given protein. Most
#'        proteins map to one uniprot id, but sometimes a protein has multiple
#'        matches.
#' @param uniprotinfo data.table with an extract of the r_uniprot table from
#'        the genoscores database. The table contains the uniprot-id, the
#'        protein name and the name of the first gene that maps to that protein
#'        name.
#' @param col.name column name whose values will be concatenated using a comma
#'        ','. This can be the gene symbol corresponding to each uniprot id, or
#'        other attributes of the gene such as chromosome and position.
#'
#' @return
#' String with concatenated values.
concat.protein <- function(ids.list, uniprotinfo, col.name) {
    ids <- unlist(strsplit(ids.list, ","))
    concat <- as.character(uniprotinfo[accessionid == ids[1], ..col.name])
    if(length(ids) > 1) {
        for(id in ids[-1]) {
            concat <- paste0(concat, ",",
                             uniprotinfo[accessionid == id, ..col.name])
        }
    }
    return(concat)
}


#' Map protein traits to their corresponding genes.
#'
#' Join tables traits, r_uniprot and r_ensembl from the genoscores database
#' to retrieve the corresponding gene symbol and gene location for all the
#' traits that are mapped to uniprot ids.
#'
#' @return
#' Data.table with gene symbol and gene location for each protein trait. If a
#' protein trait is mapped to mulitple uniprot ids, then the value in the gene
#' symbol and gene location columns is a concatenation of all corresponding
#' values separated by comma.
proteins2genes <- function() {
    require(DBI)
    scorecon <- genoscores:::getscorecon()

    ## retrieve rows from traits table mapped to a uniprot id (i.e. proteins)
    traitinfo <- dbGetQuery(scorecon, paste0(
        "select g.gwasid, g.studyid, t.* ",
        "from (select * from traits where mapped_resource = 'uniprot_id') as t ",
        "inner join gwas as g on t.traitid = g.traitid"))
    traitinfo <- as.data.table(traitinfo)

    ## retrieve mapped gene for all traits with a uniprot id (i.e. proteins)
    uniprotids <- unique(unlist(strsplit(traitinfo$mapped_value, ",")))
    uniprotinfo <- dbGetQuery(scorecon, paste0(
        "select accessionid, protein_name, first_gene_name from r_uniprot ",
        "where accessionid in ('", paste(uniprotids, collapse = "', '"),  "')"))

    ## retrieve Ensembl info for genes mapped to proteins
    geneinfo <- dbGetQuery(scorecon, paste0(
        "select * from r_ensembl where gene_symbol in ('",
        paste(unique(uniprotinfo$first_gene_name), collapse = "', '"), "')"))

    ## merge uniprot info with Ensembl info and convert to data.table
    uniprotinfo <- merge(geneinfo, uniprotinfo, by.x = "gene_symbol",
                         by.y = "first_gene_name", all = TRUE, sort = FALSE)
    uniprotinfo <- as.data.table(uniprotinfo)

    ## concatenate gene info for proteins that map to mulitple uniprot ids (and
    ## therefore multiple corresponding genes)
    traitinfo[, gene_symbol := sapply(mapped_value, function(x)
        concat.protein(x, uniprotinfo, "gene_symbol"))]
    traitinfo[, gene_chrom := sapply(mapped_value, function(x)
        concat.protein(x, uniprotinfo, "chrom"))]
    traitinfo[, gene_startpos := sapply(mapped_value, function(x)
        concat.protein(x, uniprotinfo, "startpos"))]
    traitinfo[, gene_endpos := sapply(mapped_value, function(x)
        concat.protein(x, uniprotinfo, "endpos"))]

    return(traitinfo)
}

#' Annotate scoresinfo with protein metadata
annotate.scoresinfo <- function(scoresinfo) {
    ## Get proteininfo annotation object
    proteininfo <- proteins2genes()

    idx <- scoresinfo[gwasid %in% proteininfo$gwasid, which=TRUE]
    ## index corresponding protein trait in proteininfo object
    mapidx <- match(scoresinfo$gwasid[idx], proteininfo$gwasid)
    ## fill-in scoresinfo columns for protein scores based on proteininfo
    scoresinfo[idx, gene_symbol := proteininfo$gene_symbol[mapidx]]
    scoresinfo[idx, gene_chrom := proteininfo$gene_chrom[mapidx]]
    # cols <- c("gene_startpos", "gene_endpos")
    # scoresinfo[, (cols) := lapply(.SD, as.character), .SDcols = cols]
    scoresinfo[idx, gene_startpos := proteininfo$gene_startpos[mapidx]]
    scoresinfo[idx, gene_endpos := proteininfo$gene_endpos[mapidx]]
    ## derive gene_distance for regional protein scores
    scoresinfo[gwasid %in% proteininfo$gwasid & region != 0,
               gene_distance := suppressWarnings(
                   genoscores:::distance.ranges(as.integer(gene_startpos),
                                                as.integer(gene_endpos),
                                                startpos, endpos))]
    scoresinfo[gwasid %in% proteininfo$gwasid &
               region != 0 &
               chrom_int != gene_chrom,
               gene_distance := Inf]

    ## handle gene_distance separately for proteins mapping to multiple genes
    multidx <- scoresinfo[gwasid %in% proteininfo$gwasid &
                          region != 0 &
                          grepl(",", gene_symbol), which=TRUE]

    for(jj in multidx) {
        cols <- c("gene_chrom", "gene_startpos", "gene_endpos")
        jj.info <- scoresinfo[jj, lapply(.SD, function(x) unlist(strsplit(x, ","))),
                              .SDcols = cols]
        match.chrom <- jj.info$gene_chrom == scoresinfo[jj, chrom_int]
        if(any(match.chrom)) {
            ## get minimum distance between score and gene across all mapped
            ## genes in the same chromosome as the score
            scoresinfo[jj, gene_distance := min(genoscores:::distance.ranges(
                           as.integer(jj.info$gene_startpos[match.chrom]),
                           as.integer(jj.info$gene_endpos[match.chrom]),
                           startpos, endpos))]
        } else {
            ## distance is Inf if mapped genes are not in the same chromosome
            ## as the score
            scoresinfo[jj, gene_distance := Inf]
        }
    }

    ## derive qtl_type for regional protein scores
    scoresinfo[gwasid %in% proteininfo$gwasid & region != 0,
               qtl_type := cut(gene_distance, c(0, 5e4, 5e6, Inf),
                               include=TRUE, labels=c("cis", "cis-x", "trans"))]
    scoresinfo[gwasid %in% proteininfo$gwasid & region != 0 &
                   chrom_min == 6 & startpos > 25e6 & endpos < 33e6 &
                   gene_chrom == 6 & gene_startpos > 25e6 & gene_endpos < 33e6,
               qtl_type := "cis"]

    return(scoresinfo)
}
