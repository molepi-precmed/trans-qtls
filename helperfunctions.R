select.matched <- function(x, y) {
    x <- unlist(strsplit(x, split=","))
    x <- x[x %in% y]
    x <- paste(x, collapse=",")
    return(x)
}

v.select.matched <- function(vec1, y) {
    sapply(vec1, select.matched, y)
}

#! returns a logical vector, i th element TRUE if one of the strings in vec1[i] is matched by one of the strings in vec2[i]
matchedin.elementwise <- function(vec1, vec2) { # vec1 and vec2 are vectors of comma-separated strings, of same length
    N <- length(vec1)
    matched <- logical(N)
    for(i in 1:N) {
        v1 <- unlist(strsplit(vec1[i], split=","))
        v2 <- unlist(strsplit(vec2[i], split=","))
        matched[i] <- any(v1 %in% v2)
    }
    return(matched)
}

matchedin <- function(str1, str2) { # returns TRUE if any of comma-separated strings in str1 (possibly a vector of comma-separated strings) are matched by comma-separated strings in str2
    v2 <- unlist(strsplit(str2, split=","))
    match.element <- function(str1) { # str1 is not a vector
        v1 <- unlist(strsplit(str1, split=","))
        v1.matched <- paste(unique(v1[v1 %in% v2]), collapse=",")
        return(v1.matched)
    }
    sapply(X=str1, FUN=match.element, simplify=TRUE, USE.NAMES=FALSE)
}
#matchedin(trans.genome.wide.scoresinfo[gene_symbol=="CTLA4", nearby_genes], t1d.genes$nearest)

minus.log10.phnorm.extreme <- function(z, upper=TRUE) {
    z.hi <- function(z) {
        ## https://www.johndcook.com/blog/norm-dist-bounds/
        c.upper <- 8 / pi
        c.lower <- 4
        z.hi <- z[z >= 8.2]
        x.upper <-  2 * sqrt(2 / pi) / (z.hi + sqrt(z.hi^2 + c.upper))
        x.lower <-  2 * sqrt(2 / pi) / (z.hi + sqrt(z.hi^2 + c.lower))
        ln_p.upper = -0.5 * z.hi^2 + log(x.upper)
        ln_p.lower = -0.5 * z.hi^2 + log(x.upper)
        log10p <- 0.5 * (ln_p.upper + ln_p.lower) / log(10)
        return(log10p)
    }
    z.dt <- data.table(z=abs(z), log10p=NA)
    ## phnorm gives distribution function for half-normal distribution
    z.dt[z < 8.2,  log10p := log10(1 - phnorm(z))]
    z.dt[z >= 8.2, log10p := z.hi(z)]
    return(z.dt[, -log10p])
}

p.threshold <- function(p, fdr=0.01) {
    pvalues <- data.table(p)
    pvalues[, index := .I]
    setorder(pvalues, p)
    pvalues[, i :=  .I]
    pvalues[, threshold := fdr * i / .N]
    pvalues[, is.signif := p < threshold]
    # if no pvalues are below threshold, this returns -Inf
    p.max <- suppressWarnings(signif(pvalues[is.signif==TRUE, max(p)], digits=1))
    if(is.infinite(p.max)) {
        p.max <- min(pvalues$threshold) # defaults to bonferroni threshold
    }
   return(p.max)
}

pnorm.extreme <- function(z, upper=TRUE) {
    ## https://www.johndcook.com/blog/norm-dist-bounds/
    if(upper) { # upper bound on log10 tail prob
        c <- 8 / pi
    } else { # lower bound
        c = 4
    }
    x <-  2 * sqrt(2 / pi) / (z + sqrt(z^2 + c))
    ln_p = -0.5 * z^2 + log(x)
    log10p <- ln_p / log(10)
    exponent <- floor(log10p)
    coeff <- 10^(log10p - exponent)
    string <- paste0(round(coeff), "E", exponent)
    return(string)
}

#' reformat for LaTeX p-values that are in scientific notation
format.scinot.pvalue <- function(x, nexp=1) {
    x <- toupper(x)
    x.split <- as.numeric(unlist(strsplit(as.character(x), "E"))) # split x.split at E
    x.split <- signif(as.numeric(x.split, 1))
    x.split <- t(matrix(x.split, nrow=2))
    roundedto10 <- x.split[, 1]==10
    x.split[, 1][roundedto10] <- 1
    x.split[, 2][roundedto10] <- x.split[, 2][roundedto10] - 1
    p.latex <- sprintf("\\ensuremath{%.*f \\times 10^{%0*d}}", 0, x.split[, 1], nexp, x.split[, 2])
    return(p.latex)
}

#' generate formatted p-values from z value
#' global variables: sigfig, neglogp.threshold.scinot, neglogp.threshold
format.z.aspvalue <- function(z) {
    if(!exists("sigfig")) {
        sigfig <- 1
    } else {
        if(sigfig==0) {
            sigfig <- 1
        }
    }


    ## neglogp.threshold.scinot is threshold for using scientific notation
    if(!exists("neglogp.threshold.scinot")) {
        neglogp.threshold.scinot <- 3
    } else {
        if(neglogp.threshold.scinot==0) {
            neglogp.threshold.scinot <- 3
        }
    }

    p <- signif(2 * pnorm(-abs(z)), sigfig)
    p.char <- toupper(as.character(p))
    ## pnorm.extreme returns a character string of form "NE-NNN"
    p.char[!is.na(p.char) & p.char=="0"] <- pnorm.extreme(z[!is.na(p.char) & p.char=="0"])    # where R outputs 0
    sci.revert <- grepl("E", p.char) & p > 10^-neglogp.threshold.scinot
    p.char[sci.revert] <-  format(p[sci.revert], scientific=FALSE)

    if(exists("neglogp.threshold")) {
        if(neglogp.threshold > 0) { # thresholding of p values
            p.char[p < 10^-neglogp.threshold] <- paste0("<",
                                                        format(10^-neglogp.threshold,
                                                               scientific=FALSE))
        }
    } else {
        }
    p.char[grep("E", p.char)] <- format.scinot.pvalue(p.char[grep("E", p.char)])
    return(p.char)
}

format.pvalue <- function(z, pvalue=NULL) {
    format.z.aspvalue(z)
}

## format a vector of pvalues in LaTeX and return a vector of mode character
pvalue.latex <- function(x, n=1, nexp=1, s.threshold=s.threshold, p.threshold=p.threshold) {
    ## this function has to be able to handle x whether numeric or character
    pvalue <- sapply(x, function(z) { # sapply returns a vector applying FUN to each element of x

        if (is.na(z) | is.nan(z)) {
            return(NA)
        } else if(as.numeric(z) >= 10^s.threshold) {
            ## return character string to one sig fig, not in scientific notation
            return(as.character(signif(as.numeric(z), 1)))
        } else {
            if(is.numeric(z)) {
                ## rounds to 1 sig fig and convert to character string
                z <- sprintf("%.*E", 0, signif(z, n)) # default is 1 sig fig
            } else {
                z <- toupper(z)
            }
            z <- as.numeric(unlist(strsplit(as.character(z), "E"))) # split z at E
            sprintf("\\ensuremath{%.*f\\times 10^{%0*d}}", 0, z[1], nexp, z[2])
        }
    }
    )
    pvalue <- as.character(pvalue)
    if(threshold==-4) {
        pvalue[grep("\\\\times", pvalue)] <- "<0.0001" # fix for thresholding at 0.0001
    }
    return(pvalue)
}

dist2order = function(corr, method, ...) {
  d_corr = as.dist(1 - corr)
  s = seriation::seriate(d_corr, method = method, ...)
  i = seriation::get_order(s)
  return(i)
}

diversity <- function(x, q) {
    p <- x / sum(x)
    if(q != 1) {
        diversity <- (sum(p^q))^(1 / (1 - q))
    } else {
        diversity <- 2^sum(-p * log2(p))
    }
    return(diversity)
}

minus.log10.phnorm.extreme <- function(z, upper=TRUE) {
    z.hi <- function(z) {
        ## https://www.johndcook.com/blog/norm-dist-bounds/
        c.upper <- 8 / pi
        c.lower <- 4
        z.hi <- z[z >= 8.2]
        x.upper <-  2 * sqrt(2 / pi) / (z.hi + sqrt(z.hi^2 + c.upper))
        x.lower <-  2 * sqrt(2 / pi) / (z.hi + sqrt(z.hi^2 + c.lower))
        ln_p.upper = -0.5 * z.hi^2 + log(x.upper)
        ln_p.lower = -0.5 * z.hi^2 + log(x.upper)
        log10p <- 0.5 * (ln_p.upper + ln_p.lower) / log(10)
        return(log10p)
    }
    log10p <- numeric(length(z))
    z <- abs(z)
    log10p[z < 8.2] <- log10(1 - phnorm(z[z < 8.2]))
    log10p[z >=8.2] <- z.hi(z)
    return(-log10p)
}

eqtls.gene <- function(gsymbol, qtl.info) {
    eqtls <- qtl.info[gene_symbol==gsymbol,
                                            .(qtl_type, CHR=chrom_min, startpos,
                                              endpos, variance, near.t1dgenes, nearby_genes)]
    #eqtls[, CHR := factor(CHR, levels=c(1:22, "X", "Y"))]
    setkey(eqtls, CHR)
    setkey(chr.lengths, CHR)
    eqtls <- chr.lengths[eqtls]
    eqtls[, xstart := 0.5 * (startpos + endpos) + cumbp]
    eqtls[, xend := eqtls[qtl_type=="cis", xstart]]
    return(eqtls)
}

chr.lengths <- fread(text="CHR,length
1, 248956422
2, 242193529
3,198295559
4,190214555
5,181538259
6,170805979
7,159345973
8,145138636
9,138394717
10,133797422
11,135086622
12,133275309
13,114364328
14,107043718
15,101991189
16,90338345
17,83257441
18,80373285
19,58617616
20,64444167
21,46709983
22,50818468
X,156040895
Y,57227415")
chr.lengths[, CHR :=  factor(CHR, levels=c(1:22, "X", "Y"))]
chr.lengths[, cumbp := c(0, cumsum(as.numeric(length)[-.N]))]
chr.lengths[, midpoint.chr := 0.5 * length]

manhattan <- function(results, point.size=1) {
    N <- nrow(results)
    ## thin data points
    keep.rows <- c(results[, sample(.N, ceiling(0.01 * .N))],
                   results[minuslog10p > 1, sample(.N, ceiling(0.1 * .N))],
                   results[minuslog10p > 1.5, which=TRUE])
    keep.rows <- sort(unique(keep.rows))
    results.pruned <- results[keep.rows]
    chrom.labels <- as.character(1:22)
    chrom.labels[c(19, 21)] <- " "
    chrom.midpoints <- position.absolute(chr.lengths$CHR, chr.lengths$midpoint.chr)[1:22]
    chrom.breaks <- chr.lengths[1:22, cumbp]
    chrom.ends <- position.absolute(chr.lengths$CHR, chr.lengths$length)
    p.manhattan <- ggplot(results.pruned,
                          aes(x=x, y=minuslog10p, colour=CHR)) +
        geom_point(size=point.size) +
        scale_x_continuous(expand=c(0, 0),
                           breaks=chrom.midpoints, labels=chrom.labels,
                           minor_breaks = chrom.breaks) +
        theme(panel.grid.major.x = element_blank()) +
        scale_y_continuous(expand=expansion(add = c(0, 1)),
                           breaks=c(0, 5, 10, 15, 20)) +
        labs(y = expression(paste(-log[10], " ", italic(p))),
             x = "Chromosome") +
        theme(legend.position="none") +
        theme(axis.title.x = element_text(size = 14),
              axis.text.x = element_text(size = 12),
              axis.title.y = element_text(size = 14),
              axis.text.y = element_text(size = 12)) +
        scale_color_manual(values=rep(c("#AAAAAA","#666666"), length(breaks.chr)))
    return(p.manhattan)
}

#! convert chr bp to absolute position x, using hg38 to set length of chr
position.absolute <- function(CHR, BP, padding=0) {
    positions <- data.table(CHR, BP)
    positions[, rownum := .I]
    setkey(positions, CHR)
    setkey(chr.lengths, CHR)
    positions <- chr.lengths[, .(CHR, cumbp)][positions]
    positions[, x := BP + cumbp + padding]
    setorder(positions, rownum)
    return(positions[, x])
}

## lift positions from hg19 (GRCh37) to hg38
liftover.positions <- function(CHR, BP, data.dir) {
    positions37.bed <- file.path(data.dir, "positions37.bed")
    positions38.bed <- file.path(data.dir, "positions38.bed")
    chain.file <- file.path(data.dir, "hg19ToHg38.over.chain")
    unmapped <- file.path(data.dir, "unmapped")

    positions37 <- data.table(CHR, BP=as.integer(BP))
    positions37[, rowname := sprintf("row%09d", .I)]
    positions37[, CHR := paste0("chr", CHR)]
    write.table(positions37[, .(chrom=CHR,
                                chromStart=BP-1,
                                chromEnd=BP,
                                name=rowname)],
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t",
                file=positions37.bed)
    cmd <- sprintf("liftOver %s %s %s %s",
                   positions37.bed, chain.file, positions38.bed, unmapped)
    system(cmd)
    positions38 <- fread(positions38.bed)
    positions38 <- positions38[, c(4, 2)]
    colnames(positions38) <- c("rowname", "pos.hg38")
    setkey(positions38, rowname)
    setkey(positions37, rowname)
    positions <- positions38[positions37]
    setorder(positions, rowname)
    return(positions[, pos.hg38])
}
