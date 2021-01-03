#'
#'
#'
#'
#'
#'
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
