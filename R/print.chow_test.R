#' Printing Chow tests
#'
#' This function is a method for the class \code{chow_test}.
#'
#' @param x An object of class \code{chow_test}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{chow_test}}
#'
#' @export

print.chow_test <- function(x) {
  if (!inherits(x = x, what = "chow_test")) {
    stop("The object is not of class 'chow_test'.", call. = FALSE)
  }
  cat("CHOW TESTS\n\n")
  cat("Null hypothesis: parameter constancy")
  cat("\n\nSAMPLE-SPLIT\n\n")
  cat(paste0("Chi-squared: ", formatC(x = x$stat_SS, digits = 5, format = "f"), ", df: ", x$df_SS, ", p-value: ", formatC(x = x$pval_SS, digits = 5, format = "f")))
  cat("\n\nBREAK-POINT\n\n")
  cat(paste0("Chi-squared: ", formatC(x = x$stat_BP, digits = 5, format = "f"), ", df: ", x$df_BP, ", p-value: ", formatC(x = x$pval_BP, digits = 5, format = "f")))
  cat("\n\nFORECAST\n\n")
  cat(paste0("F-value: ", formatC(x = x$stat_F, digits = 5, format = "f"), ", df: (", paste(trunc(x$df_F), collapse = ", "), "), p-value: ", formatC(x = x$pval_F, digits = 5, format = "f")))
}
