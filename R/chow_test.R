#' Chow structural break test
#'
#' This function tests the presence of a structural change by implementing the sample-split (\code{SS}), the break-point (\code{BP}) and the forecast (\code{F}) versions of the Chow test.
#'
#' @param object  An object of class \code{VECM}.
#' @param t_break An integer. The period for which the structural break is tested.
#'
#' @return An object of class \code{chow_test} is a list with the following components: \cr
#'
#' \item{stat_SS}{The sample-split test statistic.}
#' \item{df_SS}{The degrees of freedom of the test statistic \code{stat_SS}.}
#' \item{pval_SS}{The p-value of the test statistic \code{stat_SS}.}
#' \item{stat_BP}{The break-point test statistic.}
#' \item{df_BP}{The degrees of freedom of the test statistic \code{stat_BP}.}
#' \item{pval_BP}{The p-value of the test statistic \code{stat_BP}.}
#' \item{stat_F}{The forecast test statistic.}
#' \item{df_F}{The degrees of freedom of the test statistic \code{stat_F}.}
#' \item{pval_F}{The p-value of the test statistic \code{stat_F}.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Doornik, J. A. & Hendry, D. F. (1997). \emph{Modelling Dynamic Systems Using PcFiml 9.0 for Windows}, International Thomson Business Press, London.
#'
#' Lütkepohl, H., & Krätzig, M. (2004). \emph{Applied time series econometrics}. Cambridge: Cambridge University Press.
#'
#' Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(bvartools)
#'
#' ## The data for this example are taken from Lütkepohl, 2005.
#' data("e6")
#'
#' model <- VECM(data = e6, p = 3,  r = 1, method = "EGLS", spec = "uconst", season = 4)
#'
#' chow <- chow_test(object = model, t_break = 43)
#' print(chow)
#' }
#'
#' @importFrom stats pchisq pf

chow_test <- function(object, t_break) {

  if (!inherits(object, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }

  model <- object$model

  K <- model$K
  r <- model$r
  p <- model$p

  data <- model$data
  data_1 <- data[1:t_break, ]
  data_2 <- data[-(1:t_break), ]

  model_1 <- VECM(data = data_1, p = model$p, r = model$r, method = model$method, spec = model$spec, shift = model$shift, season = model$season, season_ref = model$season_ref, exogen_ect = model$exogen_ect, q_ect = model$q_ect, exogen = model$exogen, q = model$q)
  model_2 <- VECM(data = data_2, p = model$p, r = model$r, method = model$method, spec = model$spec, shift = model$shift, season = model$season, season_ref = model$season_ref, exogen_ect = model$exogen_ect, q_ect = model$q_ect, exogen = model$exogen, q = model$q)

  n <- model$n
  n_1 <- model_1$model$n
  n_2 <- model_2$model$n

  u <- object$u
  u_1 <- u[, 1:n_1]
  u_2 <- u[, (n - n_2 + 1):n]

  sigma_u <- tcrossprod(x = u, y = u) / n
  sigma_12 <- (tcrossprod(x = u_1, y = u_1) + tcrossprod(x = u_2, y = u_2)) / (n_1 + n_2)
  sigma_1_2 <- tcrossprod(x = u_1, y = u_1) / n_1 + tcrossprod(x = u_2, y = u_2) / n_2
  sigma_1 <- tcrossprod(x = model_1$u, y = model_1$u) / n_1
  sigma_2 <- tcrossprod(x = model_2$u, y = model_2$u) / n_2

  stat_SS <- (n_1 + n_2) * (log(det(sigma_12)) - log(det((n_1 * sigma_1 + n_2 * sigma_2) / (n_1 + n_2))))
  df_SS <- 2 * K * r - r^2 + p * K^2 + K
  pval_SS <- pchisq(q = stat_SS, df = df_SS, lower.tail = FALSE)

  stat_BP <- (n_1 + n_2) * log(det(sigma_1_2)) - n_1 * log(det(sigma_1)) - n_2 * log(det(sigma_2))
  df_BP <- df_SS + K * (K + 1) / 2
  pval_BP <- pchisq(q = stat_BP, df = df_BP, lower.tail = FALSE)

  R2 <- 1 - (n_1 / n)^K * det(sigma_1) / det(sigma_u)
  m <- n - n_1
  s <- sqrt(((K * m)^2 - 4) / (K^2 + m^2 - 5))
  q <- K * m / 2 - 1
  z <- n - (r + ncol(model_1$gamma)) - m - (K - m + 1) / 2
  stat_F <- (1 - (1 - R2)^(1 / s)) / (1 - R2)^(1 / s) * (z * s - q) / (K * m)
  df_1 <- K * m
  df_2 <- z * s - q
  pval_F <- pf(q = stat_F, df1 = df_1, df2 = df_2, lower.tail = FALSE)

  output <- list("stat_SS" = stat_SS, "df_SS" = df_SS, "pval_SS" = pval_SS,
                 "stat_BP" = stat_BP, "df_BP" = df_BP, "pval_BP" = pval_BP,
                 "stat_F" = stat_F, "df_F" = c(df_1, df_2), "pval_F" = pval_F)

  return(structure(.Data = output, class = "chow_test"))

}
