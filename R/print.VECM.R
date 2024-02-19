#' Printing VEC model
#'
#' This function is a method for the class \code{VECM}.
#'
#' @param x An object of class \code{VECM}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM}}
#'
#' @export

print.VECM <- function(x) {
  if (!inherits(x = x, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }
  cat("LOADING MATRIX\n\n")
  alpha <- x$alpha
  alpha <- formatC(x = alpha, digits = 8, format = "f")
  prmatrix(x = alpha, quote = FALSE, right = TRUE)
  cat("\nCOINTEGRATING MATRIX\n\n")
  beta <- x$beta
  beta <- formatC(x = beta, digits = 8, format = "f")
  prmatrix(x = beta, quote = FALSE, right = TRUE)
  cat("\nSHORT-RUN PARAMETER MATRIX\n\n")
  gamma <- x$gamma
  gamma <- formatC(x = gamma, digits = 8, format = "f")
  prmatrix(x = gamma, quote = FALSE, right = TRUE)
}
