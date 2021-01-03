#'
#'
#'
#'
#'
#'
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
