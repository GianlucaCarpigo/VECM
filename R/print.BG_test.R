#' Printing Breusch-Godfrey autocorrelation test
#'
#' This function is a method for the class \code{BG_test}.
#'
#' @param x An object of class \code{BG_test}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{BG_test}}
#'
#' @export

print.BG_test <- function(x) {
  if (!inherits(x = x, what = "BG_test")) {
    stop("The object is not of class 'BG_test'.", call. = FALSE)
  }
  cat("BREUSCH-GODFREY TEST\n\n")
  cat(paste("Null hypothesis: residuals are not autocorrelated up to lag", x$lag, "\n\n"))
  type <- x$type
  cat("Test version:", type, "\n\n")
  if (type == "LM"){
    cat(paste0("Chi-squared: ", formatC(x = x$stat, digits = 5, format = "f"), ", df: ", x$df, ", p-value: ", formatC(x = x$stat_pval, digits = 5, format = "f")))
  } else if (type == "LR") {
    cat(paste0("Chi-squared: ", formatC(x = x$stat, digits = 5, format = "f"), ", df: ", x$df, ", p-value: ", formatC(x = x$stat_pval, digits = 5, format = "f")))
    cat("\n\n")
    cat(paste0("Adjusted Chi-squared: ", formatC(x = x$stat_adj, digits = 5, format = "f"), ", df: ", x$df, ", p-value: ", formatC(x = x$stat_adj_pval, digits = 5, format = "f")))
    cat("\n---")
    cat("\nEdgeworth expansion correction")
  } else {
    cat(paste0("F-value: ", formatC(x = x$stat, digits = 5, format = "f"), ", df: (", paste(trunc(x$df), collapse = ", "), "), p-value: ", formatC(x = x$stat_pval, digits = 5, format = "f")))
  }
}
