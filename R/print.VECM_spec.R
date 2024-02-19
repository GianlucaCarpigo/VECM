#' Printing VEC model specification test
#'
#' This function is a method for the class \code{VECM_spec}.
#'
#' @param x An object of class \code{VECM_spec}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM_forecast}}
#'
#' @export

print.VECM_spec <- function(x) {
  if (!inherits(x = x, what = "VECM_spec")) {
    stop("The object is not of class 'VECM_spec'.", call. = FALSE)
  }
  cat("LIKELIHOOD-RATIO TEST FOR DETERMINISTIC SPECIFICATION\n\n")
  cat(paste0("H0: spec = ", x$spec_H0, "\nH1: spec = ", x$spec_H1, "\n\n"))
  cat(paste0("Chi-squared: ", formatC(x = x$stat, digits = 5, format = "f"), "  df: ", x$df, "  p-value: ", formatC(x = x$p_value, digits = 5, format = "f")))
}
