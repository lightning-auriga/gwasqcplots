#' Apply basic ggplot2 theme settings
#'
#' A simple utility definition of some straightforward
#' plot settings. Nothing too interesting
my.theme <- ggplot2::theme_light() + ggplot2::theme(
  plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
  axis.title = ggplot2::element_text(size = 16),
  axis.text = ggplot2::element_text(size = 14)
)
