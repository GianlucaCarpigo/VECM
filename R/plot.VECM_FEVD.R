#' Plotting VEC model FEVDs
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
#'
#' @importFrom graphics barplot legend title
#' @importFrom viridis viridis_pal

plot.VECM_FEVD <- function(x) {
  if (!inherits(x = x, what = "VECM_FEVD")) {
    stop("The object is not of class 'VECM_FEVD'.", call. = FALSE)
  }
  K <- dim(x$W)[2]
  h <- dim(x$W)[1]
  var_names <- colnames(x$W[,, 1])
  cols <- viridis_pal()(K)
  for (i in 1:K) {
    barplot(t(x$W[,, i]), names.arg = as.character(1:h), xlab = "", ylab = "", col = cols)
    title(xlab = "horizon", line = 2.4)
    title(ylab = "percentage", line = 2.4)
    title(main = paste0("variance decomposition of ", var_names[i]))
    legend(x = "top", legend = var_names, bty = "n", horiz = TRUE, xpd = TRUE, inset = c(0, 1.3), fill = cols)
  }
}
