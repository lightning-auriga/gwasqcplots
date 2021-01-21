#' Given input data and plot preferences, compute
#' color for all points and, if needed, locations for
#' labels of known loci.
#'
#' This is a helper function wrapping a bunch of
#' messy calculations about where to put known loci
#' labels and how to set point colors. This could be
#' significantly improved with integration with external
#' reference LD information, to make the locus coloring
#' anything other than a messy hack.
#'
#' @param raw.data data.frame data frame containing input
#' data for Manhattan plot
#' @param chr.col.header character vector name of chromosome
#' column header in raw.data
#' @param pos.col.header character vector name of physical
#' position column header in raw.data
#' @param rsid.col.header character vector name of rsid
#' (or other format of SNP ID) column header in raw.data
#' @param pval.col.header character vector name of p-value
#' column header in raw.data
#' @param snp.color.one character vector name of color for
#' first of alternating chromosome colors. Seems to look
#' best as one of two shades of the same color, in this case
#' grey
#' @param snp.color.two character vector name of color for
#' second of alternating chromosome colors. Seems to look
#' best as one of two shades of the same color, in this case
#' grey
#' @param prev.loci.color character vector name of color
#' with which to highlight variants (and neighbors) from
#' previously-known locus set
#' @param novel.loci.color character vector name of color
#' with which to highlight variants (and neighbors) from
#' novel locus set
#' @param filename.prev.hits character vector name of file
#' containing previously-known loci, or NA
#' @param filename.novel.hits character vector name of file
#' containing novel loci, or NA
#' @param locus.width numeric ad hoc number of variants
#' physically surrounding known variant to flag with the same
#' highlight color; this should eventually be deprecated in
#' favor of LD-based determination 
#' @return list with the following entries: snp.labels,
#' a character vector with labels to overlay on the plot denoting
#' specified loci; snp.label.x.positions, numeric vector
#' containing x-axis positions for overlay labels; snp.label.y.positions,
#' numeric vector containing y-axis positions for overlay labels;
#' all.colors, character vector containing the color assignment
#' for each point in the desired Manhattan plot
compute.manhattan.settings <- function(raw.data,
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
                                       locus.width) {
  print(paste("setting default colors for SNPs"))
  all.colors <- rep(snp.color.one, nrow(raw.data))
  all.colors[raw.data[, chr.col.header] %% 2 == 1] <- snp.color.two
  ## start the process of labelling
  snp.labels <- c()
  snp.label.x.positions <- c()
  snp.label.y.positions <- c()
  snp.label.colors <- c()
  ## for each set of previously known loci, overwrite the color
  ##   as prev.loci.color
  if (!is.na(filename.prev.hits)) {
    prev.hits <- read.table(filename.prev.hits, header = FALSE)
    prev.hits[, 1] <- as.vector(prev.hits[, 1], mode = "character")
    prev.hits[, 2] <- as.vector(prev.hits[, 2], mode = "character")
    for (i in seq_len(nrow(prev.hits))) {
      snp <- prev.hits[i, 1]
      snp.label <- prev.hits[i, 2]
      print(paste("annotating snp ", snp, sep = ""))
      if (length(raw.data[, pos.col.header][
        raw.data[, rsid.col.header] == snp
      ]) >= 1) {
        locus.pos.center <- raw.data[, pos.col.header][
          raw.data[, rsid.col.header] == snp
        ][1]
        all.colors[abs(raw.data[, pos.col.header] -
          locus.pos.center) < locus.width] <-
          prev.loci.color
        ## handle label
        snp.labels <- c(snp.labels, snp.label)
        snp.label.x.positions <- c(
          snp.label.x.positions,
          locus.pos.center
        )
        plotmax <- max(raw.data[, pval.col.header][
          abs(raw.data[, pos.col.header] - locus.pos.center) <
            locus.width
        ])
        snp.label.y.positions <- c(
          snp.label.y.positions,
          plotmax
        )
        snp.label.colors <- c(snp.label.colors, prev.loci.color)
      }
    }
  }
  ## for each set of previously unknown loci,
  ##   overwrite the color as novel.loci.color
  if (!is.na(filename.novel.hits)) {
    novel.hits <- read.table(filename.novel.hits, header = FALSE)
    novel.hits[, 1] <- as.vector(novel.hits[, 1], mode = "character")
    novel.hits[, 2] <- as.vector(novel.hits[, 2], mode = "character")
    for (i in seq_len(nrow(novel.hits))) {
      snp <- novel.hits[i, 1]
      snp.label <- novel.hits[i, 2]
      print(paste("annotating snp ", snp, sep = ""))
      stopifnot(length(raw.data[, pos.col.header][
        raw.data[, rsid.col.header] == snp
      ]) == 1)
      locus.pos.center <- raw.data[, pos.col.header][
        raw.data[, rsid.col.header] == snp
      ]
      all.colors[abs(raw.data[, pos.col.header] -
        locus.pos.center) < locus.width] <-
        novel.loci.color
      ## handle label
      snp.labels <- c(snp.labels, snp.label)
      snp.label.x.positions <- c(
        snp.label.x.positions,
        locus.pos.center
      )
      plotmax <- max(raw.data[, pval.col.header][
        abs(raw.data[, pos.col.header] - locus.pos.center) <
          locus.width
      ])
      snp.label.y.positions <- c(
        snp.label.y.positions,
        plotmax
      )
      snp.label.colors <- c(snp.label.colors, novel.loci.color)
    }
  }
  list(
    snp.labels = snp.labels,
    snp.label.x.positions = snp.label.x.positions,
    snp.label.y.positions = snp.label.y.positions,
    all.colors = all.colors
  )
}
