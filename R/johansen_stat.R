#' Johansen test statistics
#'
#' This function computes the trace and maximum eigenvalue test statistics.
#'
#' @param n_step The number of steps of a Gaussian random walk.
#' @param rw_dim The dimension of the random walk used.
#' @param spec The type of deterministic component to be considered. It can be selected between \code{none} (the default), \code{rconst}, \code{uconst}, \code{rtrend}, and \code{utrend}.
#'
#' @return The values of the trace and the maximum eigenvalue test statistics.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
#'
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats residuals

johansen_stat <- function(n_step, rw_dim, spec) {
  u <- mvrnorm(n = n_step, mu = rep(x = 0, times = rw_dim), Sigma = diag(x = 1, nrow = rw_dim, ncol = rw_dim))
  z <- matrix(data = 0, nrow = n_step, ncol = rw_dim)
  for (i in 2:n_step) {
    z[i, ] <- z[i - 1, ] + u[i - 1, ]
  }
  if (spec == "rconst") {
    z <- cbind(z, 1)
  } else if (spec == "uconst") {
    z <- residuals(lm(z ~ 1))
  } else if (spec == "rtrend") {
    z <- cbind(residuals(lm(z ~ 1)), 1:n_step - 0.5 * n_step)
  } else if (spec == "utrend") {
    x <- 1:n_step
    z <- residuals(lm(z ~ x))
  }
  D <- crossprod(x = u, y = z) %*% solve(crossprod(x = z, y = z)) %*% crossprod(x = z, y = u)
  D_trace <- sum(diag(D))
  D_max <- max(eigen(x = D, only.values = TRUE)$values)
  return(cbind(D_trace, D_max))
}
