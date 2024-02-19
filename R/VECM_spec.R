#' VEC model specification test
#'
#' This function implements a likelihood ratio test for hypotheses regarding the deterministic terms of a VEC model.
#'
#' @param data A matrix containing the multivariate time series to be analyzed.
#' @param r The number of cointegrating relations.
#' @param p The lag order of the VEC model.
#' @param spec_H0 The model with constrained deterministic terms. It can be chosen between \code{none} (the default), \code{rconst}, \code{uconst}, and \code{rtrend}.
#' @param spec_H1 The model with no constrained deterministic terms. It can be chosen between \code{rconst} (the default), \code{uconst}, \code{rtrend}, and \code{utrend}.
#' @param ... Further arguments passed to the \code{VECM} function.
#'
#' @details The model need to be tested. For example, the user cannot test \code{spec_H0 = utrend} versus \code{spec_H1 = rtrend}.
#'
#' @return An object of class \code{VECM_spec} is a list with the following components: \cr
#'
#' \item{spec_H0}{The model with constrained deterministic terms.}
#' \item{spec_H1}{The model with no constrained deterministic terms.}
#' \item{stat}{The test statistic.}
#' \item{df}{The degrees of freedom.}
#' \item{p_value}{The p-value of the test statistic.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
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
#' test_s <- VECM_spec(data = e6, r = 1, p = 3, spec_H0 = "uconst", spec_H1 = "rtrend")
#' print(test_s)
#' }
#'
#' @export

VECM_spec <- function(data, r, p, spec_H0 = c("none", "rconst", "uconst", "rtrend"), spec_H1 = c("rconst", "uconst", "rtrend", "utrend"), ...) {
  n <- nrow(data) - (p + 1)
  K <- ncol(data)
  spec_1 <- match.arg(spec_H0)
  spec_2 <- match.arg(spec_H1)
  if (spec_1 == spec_2) {
    stop("You are comparing models with the same deterministic terms.", call. = FALSE)
  }
  model_1 <- VECM(data = data, r = r, p = p, method = "ML", spec = spec_1, ...)
  model_2 <- VECM(data = data, r = r, p = p, method = "ML", spec = spec_2, ...)
  n_par_1 <- model_1$model$n_par
  n_par_2 <- model_2$model$n_par
  LR_df <- n_par_2 - n_par_1
  if (LR_df < 0) {
    stop("You are comparing models that are non-nested.", call. = FALSE)
  }
  u_1 <- model_1$u
  sigma_u_1 <- model_1$sigma_u
  log_lik_1 <- -0.5 * (n * K * log(2 * pi) + n * log(det(sigma_u_1)) + sum(diag(t(u_1) %*% solve(sigma_u_1) %*% u_1)))
  u_2 <- model_2$u
  sigma_u_2 <- model_2$sigma_u
  log_lik_2 <- -0.5 * (n * K * log(2 * pi) + n * log(det(sigma_u_2)) + sum(diag(t(u_2) %*% solve(sigma_u_2) %*% u_2)))
  LR_stat <- -2 * (log_lik_1 - log_lik_2)
  LR_pval <- pchisq(q = LR_stat, df = LR_df, lower.tail = FALSE)
  output <- list("spec_H0" = spec_1, "spec_H1" = spec_2, "stat" = LR_stat, "df" = LR_df, "p_value" = LR_pval)
  return(structure(.Data = output, class = "VECM_spec"))
}
