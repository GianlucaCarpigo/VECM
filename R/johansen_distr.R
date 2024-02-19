#' Johansen distribution
#'
#' This function simulates the Johansen distribution for the trace and the maximum eigenvalue test statistics.
#'
#' @param n_sim The number of simulations.
#' @param n_step The number of steps of a Gaussian random walk.
#' @param rw_dim The dimension of the random walk used.
#' @param spec The type of deterministic component to be considered. It can be selected between \code{none} (the default), \code{rconst}, \code{uconst}, \code{rtrend}, and \code{utrend}.
#' @param seed An optional argument specifying the seed for the simulation.
#'
#' @return A matrix containing the simulated values for the trace and maximum eigenvalue test statistics.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
#'
#' @seealso \code{\link{johansen_stat}}
#'
#' @export
#'
#' @importFrom pbapply pbreplicate
#' @importFrom MASS mvrnorm

johansen_distr <- function(n_sim, n_step, rw_dim, spec = c("none", "rconst", "uconst", "rtrend", "utrend"), seed = NULL) {
  spec <- match.arg(spec)
  if(is.null(seed)){
    set.seed(2205)
  } else {
    set.seed(seed)
  }
  distr <- t(pbreplicate(n = n_sim, expr = johansen_stat(n_step = n_step, rw_dim = rw_dim, spec = spec), simplify = "matrix"))
  colnames(distr) <- c("trace", "max_eval")
  return(distr)
}

