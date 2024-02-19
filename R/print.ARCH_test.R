#' Printing ARCH test
#'
#' This function is a method for the class \code{ARCH_test}.
#'
#' @param x An object of class \code{ARCH_test}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{ARCH_test}}
#'
#' @export

print.ARCH_test <- function(x) {
  if (!inherits(x = x, what = "ARCH_test")) {
    stop("The object is not of class 'ARCH_test'.", call. = FALSE)
  }
  cat("ARCH-LM TEST\n\n")
  cat("Null hypothesis: there is no ARCH in the residuals\n\n")
  cat("Lag:", x$lag, "\n\n")
  cat(paste0("Chi-squared: ", formatC(x = x$stat, digits = 5, format = "f"), ", df: ", x$df, ", p-value: ", formatC(x = x$stat_pval, digits = 5, format = "f")))
}

