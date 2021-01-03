#'
#'
#'
#'
#'
#'
make.qq.manhattan <- function(filename.pvals,
                                 trait.name,
                                 filename.prev.hits = NA,
                                 filename.novel.hits = NA,
                                 output.filestem = "output",
                                 gc.correct = FALSE,
                                 write.locus.labels = TRUE,
                                 rsid.header = "SNP",
                                 pval.header = "PVAL",
                                 chr.header = "CHR",
                                 pos.header = "POS",
                                 background.snp.color.one = "darkgrey",
                                 background.snp.color.two = "lightgrey",
                                 known.locus.color = "blue",
                                 novel.locus.color = "red",
                                 gws.threshold = NA,
                                 truncate.p = NA) {
  ## input error checking
  stopifnot(is.vector(filename.pvals, mode = "character"))
  stopifnot(is.vector(trait.name, mode = "character"))
  stopifnot(is.vector(filename.prev.hits, mode = "character") |
    is.na(filename.prev.hits))
  stopifnot(is.vector(filename.novel.hits, mode = "character") |
    is.na(filename.novel.hits))
  stopifnot(is.vector(output.filestem, mode = "character"))
  stopifnot(is.logical(gc.correct))
  stopifnot(is.logical(write.locus.labels))
  stopifnot(is.vector(rsid.header, mode = "character") |
    (is.na(rsid.header) &
      is.na(filename.prev.hits) &
      is.na(filename.novel.hits)))
  stopifnot(is.vector(pval.header, mode = "character"))
  stopifnot(is.vector(chr.header, mode = "character") |
    (is.na(chr.header) &
      is.na(pos.header)))
  stopifnot(is.vector(pos.header, mode = "character") |
    (is.na(chr.header) & is.na(pos.header)))
  stopifnot(is.vector(background.snp.color.one, mode = "character"))
  stopifnot(length(colours()[colours() == background.snp.color.one]) == 1)
  stopifnot(is.vector(background.snp.color.two, mode = "character"))
  stopifnot(length(colours()[colours() == background.snp.color.two]) == 1)
  stopifnot(is.vector(known.locus.color, mode = "character"))
  stopifnot(length(colours()[colours() == known.locus.color]) == 1)
  stopifnot(is.vector(novel.locus.color, mode = "character"))
  stopifnot(length(colours()[colours() == novel.locus.color]) == 1)
  ## "global" variables
  rsid.col.header <- rsid.header
  pval.col.header <- pval.header
  chr.col.header <- chr.header
  pos.col.header <- pos.header
  pos.adjust.factor <- 100
  chr.buffer <- 8000000 / pos.adjust.factor
  snp.color.one <- background.snp.color.one
  snp.color.two <- background.snp.color.two
  locus.width <- 250000 / pos.adjust.factor
  prev.loci.color <- known.locus.color
  novel.loci.color <- novel.locus.color
  qq.sim.nsims <- 1000
  qq.sim.npoints <- 10000
  ## DO NOT CHANGE THESE GLOBALS
  qq.only <- is.na(pos.header) & is.na(chr.header)
  ## read data
  print(paste("reading file ", filename.pvals, sep = ""))
  raw.data <- read.table(filename.pvals, header = TRUE)
  ## more input error checking
  stopifnot(length(colnames(raw.data)[colnames(raw.data) ==
    rsid.col.header]) == 1 |
    (is.na(filename.prev.hits) & is.na(filename.novel.hits)))
  stopifnot(length(colnames(raw.data)[colnames(raw.data) ==
    pval.col.header]) == 1)
  stopifnot(length(colnames(raw.data)[colnames(raw.data) ==
    chr.col.header]) == 1 | qq.only)
  stopifnot(length(colnames(raw.data)[colnames(raw.data) ==
    pos.col.header]) == 1 | qq.only)
  ## reduce simulations if needed
  qq.sim.npoints <- min(qq.sim.npoints, nrow(raw.data))
  ## deal with annoying factors
  raw.data[, rsid.col.header] <- as.vector(raw.data[, rsid.col.header],
    mode = "character"
  )
  if (!is.na(truncate.p)) {
    raw.data <- raw.data[!is.na(raw.data[, pval.col.header]), ]
    raw.data[, pval.col.header][raw.data[, pval.col.header] <
      truncate.p] <- truncate.p
  }
  ## determine chromosome boundaries
  if (!qq.only) {
    print(paste("recalculating position boundaries"))
    raw.data[, pos.col.header] <- raw.data[, pos.col.header] /
      pos.adjust.factor
    min.chr <- min(raw.data[, chr.col.header])
    max.chr <- max(raw.data[, chr.col.header])
    chr.lower.bound <- c()
    chr.upper.bound <- c()
    for (i in seq(min.chr, max.chr, 1)) {
      current.chr <- raw.data[, pos.col.header][
        raw.data[, chr.col.header] == i
      ]
      if (length(current.chr) > 0) {
        ## shift every chromosome so the first SNP is at position 1
        raw.data[, pos.col.header][raw.data[, chr.col.header] == i] <-
          raw.data[, pos.col.header][
            raw.data[, chr.col.header] == i
          ] - min(current.chr) + 1
        current.chr <- raw.data[, pos.col.header][
          raw.data[, chr.col.header] == i
        ]
        chr.lower.bound <- c(chr.lower.bound, min(current.chr))
        chr.upper.bound <- c(chr.upper.bound, max(current.chr))
      } else {
        chr.lower.bound <- c(chr.lower.bound, 0)
        chr.upper.bound <- c(chr.upper.bound, 0)
      }
    }
    ## adjust positions in memory to x-axis scale
    for (i in 2:length(chr.lower.bound)) {
      raw.data[, pos.col.header][raw.data[, chr.col.header] == i] <-
        raw.data[, pos.col.header][raw.data[, chr.col.header] == i] +
        chr.upper.bound[i - 1] + chr.buffer
      chr.lower.bound[i] <- chr.lower.bound[i] + chr.upper.bound[i - 1] +
        chr.buffer
      chr.upper.bound[i] <- chr.upper.bound[i] + chr.upper.bound[i - 1] +
        chr.buffer
    }
    ## determine center of each chromosome in x-axis units
    chr.center <- (chr.lower.bound + chr.upper.bound) / 2
    ## set labels
    chr.seq <- unique(raw.data[, chr.col.header])
    chr.center <- chr.center[chr.seq]
    chr.labels <- ifelse(chr.seq == 23, "X",
      ifelse(chr.seq == 24, "Y",
        ifelse(chr.seq == 25, "XY",
          ifelse(chr.seq == 26, "M", paste(chr.seq))
        )
      )
    )
  }
  ## adjust the pvalues
  raw.data[, pval.col.header] <- remap.p.values(
    raw.data[, pval.col.header],
    gc.correct
  )

  y.max <- as.double(ceiling(max(raw.data[, pval.col.header])) + 2)



  manhattan.settings <- NULL
  snp.labels <- snp.label.x.positions <- snp.label.y.positions <- NULL
  all.colors <- NULL
  if (!qq.only) {
    manhattan.settings <- compute.manhattan.settings(
      raw.data,
      chr.col.header,
      pos.col.header,
      rsid.col.header,
      pval.col.header,
      snp.color.one,
      snp.color.two,
      prev.loci.color,
      novel.loci.color,
      filename.prev.hits,
      filename.novel.hits,
      locus.width
    )
    snp.labels <- manhattan.settings[["snp.labels"]]
    snp.label.x.positions <- manhattan.settings[["snp.label.x.positions"]]
    snp.label.y.positions <- manhattan.settings[["snp.label.y.positions"]]
    all.colors <- manhattan.settings[["all.colors"]]
  }

  ## get confidence bounds for qq plot
  ## qq.sim.npoints qq.sim.nsims
  render.qq.plot(
    raw.data,
    pval.col.header,
    y.max,
    qq.sim.nsims,
    qq.sim.npoints,
    output.filestem
  )
  if (!qq.only) {
    render.manhattan.plot(
      raw.data,
      pos.col.header,
      pval.col.header,
      all.colors,
      prev.loci.color,
      novel.loci.color,
      snp.color.one,
      snp.color.two,
      gws.threshold,
      chr.center,
      chr.labels,
      y.max,
      filename.prev.hits,
      filename.novel.hits,
      write.locus.labels,
      snp.labels,
      snp.label.x.positions,
      snp.label.y.positions,
      output.filestem
    )
  }
}
