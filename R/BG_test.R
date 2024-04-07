#' Breusch-Godfrey autocorrelation test
#'
#' This function tests the serial correlation of residuals using the Breusch-Godfrey test. The user can choose between three versions of the test: Lagrange multiplier (\code{LM}), likelihood ratio (\code{LR}) and Rao.
#'
#' @param object An object of class \code{VECM}.
#' @param lag The lag order of auxiliary regression model.
#' @param type The test version. It can be selected between \code{LM} (the default), \code{LR} and \code{Rao}.
#'
#' @details This test should be applied for a \code{lag} less than the VECM lag order. For a greater lag, it is suggested to use the portmanteau test (\code{\link{port_test}}).
#'
#' @return An object of class \code{BG_test} is a list with the following components: \cr
#'
#' \item{stat}{The test statistic.}
#' \item{stat_pval}{The p-value of the test statistic.}
#' \item{stat_adj}{The adjusted test statistic (not available for \code{Rao}).}
#' \item{stat_adj_pval}{The p-value of the adjusted test statistic (not available for \code{Rao}).}
#' \item{df}{The degrees of freedom of the test statistic.}
#' \item{type}{The version of the test.}
#' \item{lag}{The lag order of auxiliary regression model.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Lütkepohl, H., & Krätzig, M. (2004). \emph{Applied time series econometrics}. Cambridge: Cambridge University Press.
#'
#' Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
#' 
#' Edgerton, D. & Shukur, G. (1999). \emph{Testing autocorrelation in a system perspective}. Econometric Reviews.
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
#' test_LR <- BG_test(object = model, lag = 2, type = "LR")
#' test_Rao <- BG_test(object = model, lag = 4, type = "Rao")
#'
#' print(test_LR)
#' print(test_Rao)
#' }
#'
#' @export
#'
#' @importFrom stats coef embed lm pchisq pf

BG_test <- function(object, lag, type = c("LM", "LR", "Rao")) {

  if (!inherits(object, what = "VECM")) {
    stop("The object is not of class 'VECM'.")
  }

  n <- object$model$n
  K <- object$model$K
  p <- object$model$p

  t_ect <- t(crossprod(x = object$beta, y = object$model$LY))
  t_LDY <- t(object$model$LDY)

  u_dip <- t(object$u)
  u_ind <- embed(x = rbind(matrix(data = 0, nrow = lag, ncol = K), u_dip), dimension = lag + 1)[, -(1:K)]

  model_H0 <- lm(u_dip ~ -1 + t_ect + t_LDY)
  
  delta <- n - nrow(coef(model_H0))
  
  e_H0 <- model_H0$residuals
  e_H1 <- lm(u_dip ~ -1 + t_ect + t_LDY + u_ind)$residuals

  sigma_H0 <- crossprod(x = e_H0, y = e_H0) / n
  sigma_H1 <- crossprod(x = e_H1, y = e_H1) / n

  type <- match.arg(type)
  
  if (type == "LM") {
    stat <- n * (K - sum(diag(sigma_H1 %*% solve(sigma_H0))))
    df <- lag * K^2
    stat_pval <- pchisq(q = stat, df = df, lower.tail = FALSE)
    output <- list("stat" = stat, "stat_pval" = stat_pval, "df" = df, "type" = type, "lag" = lag)
  }
  
  if (type == "LR") {
    stat <- n * log(det(sigma_H0) / det(sigma_H1))
    df <- lag * K^2
    stat_pval <- pchisq(q = stat, df = df, lower.tail = FALSE)
    m <- K * lag
    delta_e <- delta - m - 0.5 * (K - m + 1)
    stat_adj <- delta_e / n * stat
    stat_adj_pval <- pchisq(q = stat_adj, df = df, lower.tail = FALSE) 
    output <- list("stat" = stat, "stat_pval" = stat_pval, "stat_adj" = stat_adj, "stat_adj_pval" = stat_adj_pval, "df" = df, "type" = type, "lag" = lag)
  }

  if (type == "Rao") {
    m <- K * lag
    s <- sqrt((K^2 * m^2 - 4) / (K^2 + m^2 - 5))
    q <- 0.5 * K * m - 1
    delta_e <- delta - m - 0.5 * (K - m + 1)
    stat <- ((det(sigma_H0) / det(sigma_H1))^(1 / s) - 1) * (delta_e * s - q) / (K * m)
    df1 <- lag * K^2
    df2 <- delta_e * s - q
    stat_pval <- pf(q = stat, df1 = df1, df2 = df2, lower.tail = FALSE)
    output <- list("stat" = stat, "stat_pval" = stat_pval, "df" = c(df1, df2), "type" = type, "lag" = lag)
  }

  return(structure(.Data = output, class = "BG_test"))
}
