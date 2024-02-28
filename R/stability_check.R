#' Stability condition check
#'
#' This function checks the eigenvalues stability condition in a VEC model.
#'
#' @param object An object of class \code{VECM}.
#' @param plot A logical value. If \code{TRUE} (the default), a plot shows the roots of the companion matrix.
#'
#' @return A list containing the following components: \cr
#'
#' \item{Re}{The real part of the roots of the companion matrix.}
#' \item{Im}{The complex part of the roots of the companion matrix.}
#' \item{Mod}{The modulus of the roots of the companion matrix.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Lütkepohl, H., & Krätzig, M. (2004). \emph{Applied time series econometrics}. Cambridge: Cambridge University Press.
#'
#' Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
#'
#' @seealso \code{\link{VECM}}
#'
#' @examples
#' \dontrun{
#' require(bvartools)
#'
#' ## The data for this example are taken from Lütkepohl, 2005.
#' data("e6")
#'
#' model <- VECM(data = e6, p = 3,  r = 1, method = "EGLS", spec = "uconst", season = 4)
#' check <- stability_check(model)
#' }
#'
#' @export
#'
#' @importFrom graphics abline lines par

stability_check <- function(object, plot = TRUE) {

  var <- VECM_to_VAR(object)

  A <- var$A
  K <- nrow(A)
  p <- ncol(A)/K
  id <- cbind(diag(x = 1, nrow = K * (p - 1), ncol = K * (p - 1)), matrix(data = 0, nrow = K * (p - 1), ncol = K))
  roots <- as.complex(eigen(rbind(A, id))$values)

  if (plot) {
    par(pty = "s")
    plot(roots, ylim = c(-1, 1), xlim = c(-1, 1), asp = 1, ylab = "Im(eigen)", xlab = "Re(eigen)", main = "Eigenvalues of the companion matrix",
    panel.first = c(lines(complex(modulus = 1, argument = 0.01 * 2 * pi)^(0:100), col = 'grey'),
    abline(h = 0, col = 'grey'), abline(v = 0, col = 'grey')))
    par(pty = "m")
  }

  return(list("Re" = Re(roots), "Im" = Im(roots), "Mod" = Mod(roots)))
}
