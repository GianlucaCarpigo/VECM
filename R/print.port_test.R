#' Printing portmanteau autocorrelation test
#'
#' This function is a method for the class \code{port_test}.
#'
#' @param x An object of class \code{port_test}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{port_test}}
#'
#' @export

print.port_test <- function(x) {
  if (!inherits(x, what = "port_test")) {
    stop("The object is not of class 'VECM_porttest'.")
  }
  cat("PORTMANTEAU TEST\n\n")
  cat(paste("Null hypothesis: residuals are not autocorrelated up to lag", x$lag, "\n\n"))
  cat(paste0("Chi-squared: ", formatC(x = x$stat, digits = 5, format = "f"), ", df: ", x$df, ", p-value: ", formatC(x = x$stat_pval, digits = 5, format = "f")))
  cat("\n\n")
  cat(paste0("Adjusted Chi-squared: ", formatC(x = x$stat_adj, digits = 5, format = "f"), ", df: ", x$df, ", p-value: ", formatC(x = x$stat_adj_pval, digits = 5, format = "f")))
  invisible(x)
}
