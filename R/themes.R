#' Apply basic ggplot2 theme settings
#'
#' A simple utility definition of some straightforward
#' plot settings. Nothing too interesting
my.theme <- theme_light() + theme(
  plot.title = element_text(size = 18, hjust = 0.5),
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14)
)
