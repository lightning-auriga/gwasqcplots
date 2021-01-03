#'
#'
#'
#'
#'
#'
render.qq.plot <- function(raw.data,
                           pval.col.header,
                           y.max,
                           qq.sim.nsims,
                           qq.sim.npoints,
                           output.filestem) {
  print(paste("calculating confidence interval for qq plot"))
  qq.sim.data <- c()
  for (i in seq_len(qq.sim.nsims)) {
    qq.sim.curdata <- sort(-log10(runif(qq.sim.npoints, 0, 1) *
      qq.sim.npoints / nrow(raw.data)))
    if (is.vector(qq.sim.data)) {
      qq.sim.data <- rbind(qq.sim.curdata)
    } else {
      qq.sim.data <- rbind(qq.sim.data, qq.sim.curdata)
    }
  }
  qq.poly.x <- sort(-log10(seq(1, qq.sim.npoints, 1) / nrow(raw.data)))
  qq.poly.x <- c(qq.poly.x, -log10(1 / nrow(raw.data)))
  qq.poly.upper.y <- c()
  qq.poly.lower.y <- c()
  for (i in seq_len(ncol(qq.sim.data))) {
    qq.poly.upper.y <- c(
      qq.poly.upper.y,
      sort(qq.sim.data[, i])[
        ceiling(qq.sim.nsims * 0.975)
      ]
    )
    qq.poly.lower.y <- c(
      qq.poly.lower.y,
      sort(qq.sim.data[, i])[
        ceiling(qq.sim.nsims * 0.025)
      ]
    )
  }
  qq.poly.upper.y <- c(qq.poly.upper.y, -log10(1 / nrow(raw.data)))
  qq.poly.lower.y <- c(qq.poly.lower.y, -log10(1 / nrow(raw.data)))
  qq.y.max <- max(c(
    y.max - 1, ceiling(max(qq.poly.upper.y)),
    ceiling(max(qq.poly.lower.y))
  ))

  ## plot the sucker
  print(paste("rendering qq plot"))
  qq.limits <- c(0, max(ceiling(-log10(1 / nrow(raw.data))), qq.y.max))
  qq.labels <- c()
  for (i in 0:qq.limits[2]) {
    qq.labels <- c(qq.labels, paste(i))
  }
  qq.data <- data.frame(
    x = sort(-log10(seq(1, nrow(raw.data), 1) / nrow(raw.data))),
    y = sort(raw.data[, pval.col.header])
  )
  qq.poly <- data.frame(
    x = c(qq.poly.x, qq.poly.x[rev(seq_along(qq.poly.x))]),
    y = c(qq.poly.lower.y, qq.poly.upper.y[rev(seq_along(qq.poly.upper.y))])
  )
  my.plot <- ggplot2::ggplot(ggplot2::aes(
    x = .data$x,
    y = .data$y
  ),
  data = qq.data
  )
  my.plot <- my.plot + my.theme + ggplot2::scale_x_continuous(
    breaks = 0:qq.limits[2],
    labels = qq.labels
  ) +
    ggplot2::scale_y_continuous(breaks = 0:qq.limits[2], labels = qq.labels)
  my.plot <- my.plot + ggplot2::geom_polygon(ggplot2::aes(
    x = .data$x,
    y = .data$y
  ),
  data = qq.poly,
  fill = "purple", alpha = 0.2
  )
  my.plot <- my.plot + ggplot2::geom_abline(slope = 1, intercept = 0)
  my.plot <- my.plot + ggplot2::geom_point()
  my.plot <- my.plot + ggplot2::xlab(expression("-log"[10] * "(Expected P)")) +
    ggplot2::ylab(expression("-log"[10] * "(Observed P)"))

  options(bitmapType = "cairo", device = "jpeg")
  ggsave(paste(output.filestem, ".qq.jpg", sep = ""),
    plot = my.plot,
    height = 10,
    width = 16 / 9 * 10,
    units = "in"
  )
}
