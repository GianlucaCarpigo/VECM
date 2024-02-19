#' Printing VEC model FEVDs
#'
#' This function is a method for the class \code{VECM_FEVD}.
#'
#' @param x An object of class \code{VECM_FEVD}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM_FEVD}}
#'
#' @export

print.VECM_FEVD <- function(x) {
  if (!inherits(x = x, what = "VECM_FEVD")) {
    stop("The object is not of class 'VECM_FEVD'.", call. = FALSE)
  }
  K <- dim(x$W)[2]
  h <- dim(x$W)[1]
  var_names <- colnames(x$W[,, 1])
  for(i in 1:K){
    cat("Variance decomposition of", var_names[i], "\n\n")
    mat <- cbind(1:h, x$forecast_error[, i], x$W[,, i])
    prmatrix(x = mat, rowlab = rep(x = "", length = h), collab = c("horizon", "error", var_names))
    if (i != K) {
      cat("\n")
    }
  }
}
