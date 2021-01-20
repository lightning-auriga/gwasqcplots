#' Convert raw p-values to plottable values with optionally
#' correction.
#'
#' Take a vector of p-values from an input file.
#' Apply GC correction if requested, handling small
#' p-values as needed. Return -log10 of either corrected
#' or uncorrected values for plotting in QQ or Manhattan plots.
#'
#' GC (genomic control) correction is an old ad hoc method
#' of adjusting inflation in the test statistic distribution
#' of a GWAS. The method converts the p-values to corresponding
#' 1-df Chi^2 test statistics, and then adjusts them all by
#' a normalization factor (median(observed) /
#' median(uniform random p-values)).
#'
#' Note that a GC inflation factor less than 1, indicating
#' test statistic *deflation*, is very possible for various,
#' usually bad, reasons, and conventionally such a correction
#' factor is unconditionally set to 1, to avoid artificially
#' inflating a distribution.
#'
#' A number of assumptions underpin the sane usage of GC
#' correction. Note most importantly that:
#'   - inflation does not necessarily mean stratification,
#'     nor does it necessarily mean polygenicity
#'   - deflation often means convergence problems with model
#'   - LD pruning should be applied before inflation factor
#'     calculation; this is a candidate for future features
#'   - it can be difficult to visually evaluate the inflation
#'     factor's suitability in a distribution due to how heavily
#'     the median is weighted by the bulk of the distribution
#'   - inflation factor calculation assumes genome-wide variants
#'     are present, among other things. do not apply GC correction
#'     when only a partial set of variants is available
#'
#' @param input.p numeric vector of input p-values
#' @param gc.correct logical; whether to apply GC correction
#' @return numeric vector of p-values, optionally GC corrected,
#' with -log10 transform applied
remap.p.values <- function(input.p,
                           gc.correct) {
  y <- input.p
  if (gc.correct) {
    print(paste("gc correcting"))
    y.low <- y[y < 1e-10]
    y.hi <- y[y >= 1e-10]
    if (length(y.low) > 0) {
      y.low <- qchisq(log(y.low), 1, lower.tail = FALSE, log.p = TRUE)
    }
    if (length(y.hi) > 0) {
      y.hi <- qchisq(y.hi, 1, lower.tail = FALSE)
    }
    inflate <- median(c(y.hi, y.low)) / 0.455
    print(paste("raw inflation factor is ", inflate, sep = ""))
    if (inflate < 1) {
      print(paste("WARNING: inflation factor less than 1 detected, ",
        "setting to 1: ", inflate,
        sep = ""
      ))
      inflate <- 1
    }
    if (length(y.low) > 0) {
      y[input.p < 1e-10] <- exp(pchisq(y.low / inflate,
        1,
        lower.tail = FALSE,
        log.p = TRUE
      ))
    }
    if (length(y.hi) > 0) {
      y[input.p >= 1e-10] <- pchisq(y.hi / inflate,
        1,
        lower.tail = FALSE
      )
    }
  }
  print(paste("remapping pvalues"))
  -log10(y)
}
