#' Multivariate normality tests
#'
#' This function performs a series of normality tests for the residuals of a VEC model.
#'
#' @param object An object of class \code{VECM}.
#' @param type The factorization of the residuals covariance matrix. It can be selected between \code{Lutkepohl} (the default), \code{Doornik-Hansen}, and \code{Urzua}.
#'
#' @return An object of class \code{VECM} is a list with the following components: \cr
#'
#' \item{b1}{The estimated skewnesses of the residuals.}
#' \item{b1_chisq}{The skewness test statsitics.}
#' \item{b1_pval}{The p-value of the test statistic \code{b1_chisq}.}
#' \item{b1_joint}{The joint skewness test statsitic.}
#' \item{b1_joint_pval}{The p-value of the test statistic \code{b1_joint}.}
#' \item{b2}{The estimated kurtosis of the residuals.}
#' \item{b2_chisq}{The kurtosis test statsitic.}
#' \item{b2_pval}{The p-value of the test statistic \code{b2_chisq}.}
#' \item{b2_joint}{The joint kurtosis test statsitic.}
#' \item{b2_joint_pval}{The p-value of the test statistic \code{b2_joint}.}
#' \item{JB}{The Jarque-Bera test statsitic.}
#' \item{JB_pval}{The p-value of the test statistic \code{JB}.}
#' \item{JB_joint}{The joint Jarque-Bera test statsitic.}
#' \item{JB_joint_pval}{The p-value of the test statistic \code{JB_joint}.}
#' \item{type}{The factorization of the residuals covariance matrix.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@uniroma1.it}
#'
#' @references Urzua, C. M. (1997). \emph{Omnibus tests for multivariate normality based on a class of maximum entropy distributions}. In Applying Maximum Entropy to Econometric Problems. Emerald Group Publishing Limited.
#'
#' Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
#'
#' Doornik, J. A., & Hansen, H. (2008). \emph{An omnibus test for univariate and multivariate normality}. Oxford bulletin of economics and statistics, 70, 927-939.
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
#' test <- normality_test(object = model, type = "Doornik-Hansen")
#'
#' print(test)
#' }
#'
#' @export
#'
#' @importFrom expm sqrtm
#' @importFrom stats cov2cor pchisq

normality_test <- function(object, type = c("Lutkepohl", "Doornik-Hansen", "Urzua")) {

  if (!inherits(object, what = "VECM")) {
    stop("The object is not of class 'VECM'.")
  }

  K <- object$model$K
  n <- object$model$n
  u <- object$u

  u_cent <- u - rowMeans(u)
  sigma_u_cent <- tcrossprod(x = u_cent, y = u_cent) / n

  type <- match.arg(type)

  if (type == "Lutkepohl") {
    u_ort <- solve(t(chol(sigma_u_cent))) %*% u_cent

    b1 <- rowMeans(u_ort^3)
    b1_chisq <- n / 6 * b1^2
    b1_pval <- pchisq(q = b1_chisq, df = 1, lower.tail = FALSE)
    b1_joint <- sum(b1_chisq)
    b1_joint_pval <- pchisq(q = b1_joint, df = K, lower.tail = FALSE)

    b2 <- rowMeans(u_ort^4)
    b2_chisq <- n / 24 * (b2 - 3)^2
    b2_pval <- pchisq(q = b2_chisq, df = 1, lower.tail = FALSE)
    b2_joint <- sum(b2_chisq)
    b2_joint_pval <- pchisq(q = b2_joint, df = K, lower.tail = FALSE)

    JB <- b1_chisq + b2_chisq
    JB_pval <- pchisq(q = JB, df = 2, lower.tail = FALSE)
    JB_joint <- sum(JB)
    JB_joint_pval <- pchisq(q = JB_joint, df = 2 * K, lower.tail = FALSE)
  }

  if (type == "Doornik-Hansen") {
    R_dec <- eigen(cov2cor(sigma_u_cent))
    eval <- R_dec$values
    evec <- R_dec$vectors
    dimnames(evec) <- dimnames(sigma_u_cent)

    u_ort <- evec %*% diag(x = eval^(-0.5), nrow = K, ncol = K) %*% t(evec) %*% diag(x = diag(sigma_u_cent)^(-0.5), nrow = K, ncol = K) %*% u_cent

    b1 <- rowMeans(u_ort^3)
    b2 <- rowMeans(u_ort^4)

    beta <- 3 * (n + 1) * (n + 3) * (n^2 + 27 * n - 70) / ((n - 2) * (n + 5) * (n + 7) * (n + 9))
    w <- sqrt(2 * (beta - 1)) - 1
    delta <- 1 / sqrt(log(sqrt(w)))
    y <- b1 * sqrt((w - 1) * (n + 1) * (n + 3) / (12 * (n - 2)))
    z1 <- delta * log(y + sqrt(y^2 + 1))
    b1_chisq <- z1^2
    b1_pval <- pchisq(q = b1_chisq, df = 1, lower.tail = FALSE)
    b1_joint <- sum(b1_chisq)
    b1_joint_pval <- pchisq(q = b1_joint, df = K, lower.tail = FALSE)

    delta <- (n - 3) * (n + 1) * (n^2 + 15 * n - 4)
    a <- (n - 2) * (n + 5) * (n + 7) * (n^2 + 27 * n - 70) / (6 * delta)
    c <- (n^2 - 49) * (n + 5) * (n^2 + 2 * n - 5) / (6 * delta)
    k <- (n + 5) * (n + 7) * (n^3 + 37 * n^2 + 11 * n - 313) / (12 * delta)
    alpha <- a + b1^2 * c
    chi <- (b2 - b1^2 - 1) * 2 * k
    z2 <- ((chi / (2 * alpha))^(1 / 3) + 1 / (9 * alpha) - 1) * sqrt(9 * alpha)
    b2_chisq <- z2^2
    b2_pval <- pchisq(q = b2_chisq, df = 1, lower.tail = FALSE)
    b2_joint <- sum(b2_chisq)
    b2_joint_pval <- pchisq(q = b2_joint, df = K, lower.tail = FALSE)

    JB <- b1_chisq + b2_chisq
    JB_pval <- pchisq(q = JB, df = 2, lower.tail = FALSE)
    JB_joint <- sum(JB)
    JB_joint_pval <- pchisq(q = JB_joint, df = 2 * K, lower.tail = FALSE)
  }

  if (type == "Urzua") {
    u_ort <- solve(sqrtm(sigma_u_cent)) %*% u_cent
    rownames(u_ort) <- rownames(sigma_u_cent)

    b1 <- rowMeans(u_ort^3)
    b2 <- rowMeans(u_ort^4)

    V_b1 <- 6 * (n - 2) / ((n + 1) * (n + 3))
    E_b2 <- 3 * (n - 1) / (n + 1)
    V_b2 <- 24 * n * (n - 2) * (n - 3) / ((n + 1)^2 * (n + 3) * (n + 5))

    b1_chisq <- b1^2 / V_b1
    b1_pval <- pchisq(q = b1_chisq, df = 1, lower.tail = FALSE)
    b1_joint <- sum(b1_chisq)
    b1_joint_pval <- pchisq(q = b1_joint, df = K, lower.tail = FALSE)

    b2_chisq <- (b2 - E_b2)^2 / V_b2
    b2_pval <- pchisq(q = b2_chisq, df = 1, lower.tail = FALSE)
    b2_joint <- sum(b2_chisq)
    b2_joint_pval <- pchisq(q = b2_joint, df = K, lower.tail = FALSE)

    JB <- b1_chisq + b2_chisq
    JB_pval <- pchisq(q = JB, df = 2, lower.tail = FALSE)
    JB_joint <- sum(JB)
    JB_joint_pval <- pchisq(q = JB_joint, df = 2 * K, lower.tail = FALSE)
  }

  output <- list("b1" = b1, "b1_chisq" = b1_chisq, "b1_pval" = b1_pval, "b1_joint" = b1_joint, "b1_joint_pval" = b1_joint_pval,
                 "b2" = b2, "b2_chisq" = b2_chisq, "b2_pval" = b2_pval, "b2_joint" = b2_joint, "b2_joint_pval" = b2_joint_pval,
                 "JB" = JB, "JB_pval" = JB_pval, "JB_joint" = JB_joint, "JB_joint_pval" = JB_joint_pval, "type" = type)

  return(structure(.Data = output, class = "normality_test"))

}
