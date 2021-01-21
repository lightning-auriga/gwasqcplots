#' Given input data and a variety of settings,
#' create a Manhattan plot and save it to file.
#'
#' This function only handles direct interactions
#' with `ggplot2`; it doesn't do much in the way of
#' input data processing. most of that is handled
#' by other library functions. the description
#' of this process will be expanded in later versions.
#'
#' @param raw.data data.frame; input data for
#' Manhattan plot with appropriate transformations
#' @param pos.col.header character vector; header
#' name for column of transformed physical positions
#' in input data
#' @param pval.col.header character vector; header
#' name for column of transformed p-values in input
#' data
#' @param all.colors character vector; color settings
#' for each row of input data
#' @param prev.loci.color character vector; color
#' setting for known loci highlighting
#' @param novel.loci.color character vector; color
#' setting for novel loci highlighting
#' @param snp.color.one character vector;
#' darker color for emphasizing alternating
#' chromosomes in Manhattan plot
#' @param snp.color.two character vector;
#' lighter color for emphasizing alternating
#' chromosomes in Manhattan plot
#' gws.threshold numeric; p-value threshold
#' to use for genome-wide significance threshold;
#' defaults to 5e-8. if no threshold should be rendered,
#' set to NA
#' chr.center numeric vector; location of chromosome
#' centers along the x-axis for labeling purposes
#' chr.labels character vector; labels to render along
#' x-axis (e.g.: "1", "2", "23", "XY" if you want, etc.)
#' y.max numeric; override value to set as maximum yintercept
#' position in output plot. should exceed valid space 
#' of observed points by some amount for aesthetically
#' pleasing plotting 
#' filename.prev.hits character vector; name of file
#' of known signals, or NA. this is only used to control
#' whether other data are processed as if these signals
#' exist; the file itself is not read. this is a candidate
#' for removal from the function interface
#' filename.novel.hits character vector; name of file
#' of novel signals, or NA. this is only used to control
#' whether other data are processed as if these signals
#' exist; the file itself is not read. this is a candidate
#' for removal from the function interface
#' write.locus.labels logical; whether to
#' render labels for known/novel signal highlights
#' snp.labels character vector; labels to be applied
#' to highlighted loci in plot. probably should be
#' renamed to "loci.labels" for consistency
#' snp.label.x.positions numeric vector; x coordinates
#' for loci labels. should probably be renamed to
#' "loci.label.x.positions" for consistency
#' snp.label.y.positions numeric vector; y coordinates
#' for loci labels. should probably be renamed to
#' "loci.label.y.positions" for consistency
#' output.filestem character vector; prefix for output
#' plot, to which will be appended ".manhattan.jpg"
render.manhattan.plot <- function(raw.data,
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
                                  output.filestem) {
    print(paste("rendering manhattan plot"))
  manhattan.data <- data.frame(
    x = raw.data[, pos.col.header],
    y = raw.data[, pval.col.header],
    colour = all.colors,
    colour.factor = factor(ifelse(all.colors == prev.loci.color,
      "Known Locus",
      ifelse(all.colors == novel.loci.color,
        "Novel Locus",
        ifelse(all.colors == snp.color.one,
          "Null Variant 1", "Null Variant 2"
        )
      )
    ))
  )
  gg.colour.levels <- c(snp.color.one, snp.color.two)
  if (length(which(all.colors == novel.loci.color)) > 0) {
    gg.colour.levels <- c(novel.loci.color, gg.colour.levels)
  }
  if (length(which(all.colors == prev.loci.color)) > 0) {
    gg.colour.levels <- c(prev.loci.color, gg.colour.levels)
  }
  my.plot <- ggplot2::ggplot(ggplot2::aes(
    x = .data$x,
    y = .data$y,
    colour = colour.factor
  ),
  data = manhattan.data
  )
  my.plot <- my.plot + my.theme +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 11),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y.right = ggplot2::element_blank()
    )
  my.plot <- my.plot + ggplot2::geom_point()
  if (!is.na(gws.threshold)) {
    my.plot <- my.plot +
      ggplot2::geom_hline(
        yintercept = -log10(gws.threshold),
        colour = "red",
        lty = 2,
        alpha = 0.2
      )
  }
  if (length(gg.colour.levels < 3)) {
    my.plot <- my.plot +
      ggplot2::scale_colour_manual(values = gg.colour.levels, guide = FALSE)
  } else {
    legend.breaks <- c()
    if (length(which(colour.factor == "Known Locus")) > 0) {
      legend.breaks <- c(legend.breaks, "Known Locus")
    }
    if (length(which(colour.factor == "Novel Locus")) > 0) {
      legend.breaks <- c(legend.breaks, "Novel Locus")
    }
    my.plot <- my.plot + ggplot2::scale_colour_manual(
      breaks = legend.breaks,
      values = gg.colour.levels,
      name = "Locus Types"
    )
  }
  my.plot <- my.plot + ggplot2::xlab("Chromosomes") +
    ggplot2::ylab(expression("-log"[10] * "(P)"))
  my.plot <- my.plot + ggplot2::scale_x_continuous(
    breaks = chr.center,
    labels = chr.labels
  )
  my.plot <- my.plot + ggplot2::scale_y_continuous(breaks = 0:y.max)

    if (!is.na(filename.novel.hits)) {
    novel.hit.data <- data.frame(
      x = raw.data[, pos.col.header][all.colors == novel.loci.color],
      y = raw.data[, pval.col.header][all.colors == novel.loci.color]
    )
    my.plot <- my.plot + ggplot2::geom_point(ggplot2::aes(
      x = .data$x,
      y = .data$y
    ),
    data = novel.hit.data,
    colour = novel.loci.color
    )
  }
    if (!is.na(filename.prev.hits)) {
    prev.hit.data <- data.frame(
      x = raw.data[, pos.col.header][all.colors == prev.loci.color],
      y = raw.data[, pval.col.header][all.colors == prev.loci.color]
    )
    my.plot <- my.plot + ggplot2::geom_point(ggplot2::aes(
      x = .data$x,
      y = .data$y
    ),
    data = prev.hit.data,
    colour = prev.loci.color
    )
  }
  if (write.locus.labels) {
    annotate("text",
      snp.label.x.positions,
      snp.label.y.positions,
      label = snp.labels
    )
  }
  options(bitmapType = "cairo", device = "jpeg")
  ggsave(paste(output.filestem, ".manhattan.jpg", sep = ""),
    plot = my.plot, height = 10, width = 16 / 9 * 10, units = "in"
    )
}
