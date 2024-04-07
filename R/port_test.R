#' Portmanteau autocorrelation test
#'
#' This function tests the serial correlation of residuals using a portmanteau test.
#'
#' @param object An object of class \code{VECM}.
#' @param lag The number of tested lagged autocorrelations.
#'
#' @return An object of class \code{port_test} is a list with the following components: \cr
#'
#' \item{stat}{The test statistic.}
#' \item{stat_pval}{The p-value of the test statistic.}
#' \item{stat_adj}{The adjusted test statistic.}
#' \item{stat_adj_pval}{The p-value of the adjusted test statistic.}
#' \item{lag}{The number of tested lagged autocorrelations.}
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
#' test_p <- port_test(object = model, lag = 10)
#'
#' print(test_p)
#' }
#'
#' @export
#'
#' @importFrom stats pchisq

port_test <- function(object, lag = 10) {

  if (!inherits(object, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }
  
  if (!is.null(object$model$exogen) | !is.null(object$model$exogen_ect)) {
    stop("The portmanteau test can not be applied with exogenous variables.", call. = FALSE)
  }

  K <- object$model$K
  N <- object$model$N
  n <- object$model$n
  p <- object$model$p
  r <- object$model$r
  
  if (lag <= p) {
    stop("The number of tested lagged autocorrelations must be greater than 'p'.", call. = FALSE)
  }
  
  u <- object$u
  u_cent <- u - rowMeans(u)

  C <- array(data = 0, dim = c(K, K, lag))

  for (i in 1:lag) {
    for (j in (i + 1):n) {
      C[,, i] <- C[,, i] + tcrossprod(x = u_cent[, j], y = u_cent[, j - i])
    }
  }

  C_lag <- matrix(data = C, nrow = lag * K^2, ncol = 1)
  C_0 <- tcrossprod(x = u_cent, y = u_cent)

  stat <- n * t(C_lag) %*% solve(diag(x = 1, nrow = lag, ncol = lag) %x% C_0 %x% C_0) %*% C_lag
  stat_adj <- n^2 * t(C_lag) %*% solve(diag(x = (n - 1):(n - lag), nrow = lag, ncol = lag) %x% C_0 %x% C_0) %*% C_lag

  df <- K^2 * (lag - p) - r * K
  stat_pval <- pchisq(q = stat, df = df, lower.tail = FALSE)
  stat_adj_pval <- pchisq(q = stat_adj, df = df, lower.tail = FALSE)

  output <- list("stat" = stat, "stat_pval" = stat_pval, "stat_adj" = stat_adj, "stat_adj_pval" = stat_adj_pval, "lag" = lag, "df" = df)

  return(structure(.Data = output, class = "port_test"))

}
