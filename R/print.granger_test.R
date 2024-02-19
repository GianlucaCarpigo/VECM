#' Printing Granger-causality test
#'
#' This function is a method for the class \code{granger_test}.
#'
#' @param x An object of class \code{granger_test}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{granger_test}}
#'
#' @export

print.granger_test <- function(x) {
  if (!inherits(x = x, what = "granger_test")) {
    stop("The object is not of class 'granger_test'.", call. = FALSE)
  }
  dip <- x$dip
  indip <- x$indip
  cat("GRANGER-CAUSALITY TEST\n\n")
  if (length(dip) == 1) {
    if (length(indip) == 1) {
      cat(paste0("Null hypothesis: ", indip, " does not Granger-cause ", dip))
    } else {
      cat(paste0("Null hypothesis: (", paste(indip, collapse = ", "), ") do not Granger-cause ", dip))
    }
  } else {
    if (length(indip) == 1) {
      cat(paste0("Null hypothesis: ", indip, " does not Granger-cause (", paste(dip, collapse = ", "), ")"))
    } else {
      cat(paste0("Null hypothesis: (", paste(indip, collapse = ", "), ") do not Granger-cause (", paste(dip, collapse = ", "), ")"))
    }
  }
  cat("\n\n")
  cat(paste0("Chi-squared: ", formatC(x = x$stat, digits = 5, format = "f"), ", df: ", x$df, ", p-value: ", formatC(x = x$stat_pval, digits = 5, format = "f")))
}
