#' ARCH test
#'
#' This function tests if ARCH effects are present in the residuals using a Lagrange multiplier test.
#'
#' @param object An object of class \code{VECM}.
#' @param lag The lag order of auxiliary regression model.
#'
#' @return An object of class \code{ARCH_test} is a list with the following components: \cr
#'
#' \item{lag}{The lag order of auxiliary regression model.}
#' \item{stat}{The test statistic.}
#' \item{stat_pval}{The p-value of the test statistic.}
#' \item{df}{The degrees of freedom of the test statistic.}
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
#' test <- ARCH_test(object = model, lag = 2)
#' print(test)
#' }
#'
#' @export
#'
#' @importFrom stats embed lm pchisq

ARCH_test <- function(object, lag = 5) {

  if (!inherits(x = object, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }

  K <- object$model$K
  N <- object$model$N

  n <- object$model$n

  u <- t(scale(t(object$u)))

  B_dim <- 0.5 * K * (K + 1)

  u_mat <- matrix(data = NA, nrow = N, ncol = B_dim)

  for (i in seq_len(n)) {
    temp <- tcrossprod(x = u[, i], y = u[, i])
    u_mat[i, ] <- temp[lower.tri(x = temp, diag = TRUE)]
  }

  Z <- embed(x = u_mat, dimension = lag + 1)

  temp <- seq.int(from = 1, to = B_dim)

  u_dip <- Z[, temp]

  e_H0 <- lm(u_dip ~ 1)$residuals
  sigma_H0 <- crossprod(x = e_H0, y = e_H0) / n

  u_ind <- Z[, -temp]

  e_H1 <- lm(u_dip ~ u_ind)$residuals
  sigma_H1 <- crossprod(x = e_H1, y = e_H1) / n

  stat <- nrow(u_dip) * B_dim * (1 - sum(diag(sigma_H1 %*% solve(sigma_H0))) / B_dim)
  df <- lag * B_dim^2
  pval <- pchisq(q = stat, df = df, lower.tail = FALSE)

  output <- list("lag" = lag, "stat" = stat, "stat_pval" = pval, "df" = df)

  return(structure(.Data = output, class = "ARCH_test"))

}


