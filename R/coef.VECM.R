#' Extract VECM coefficients
#'
#' This function is a method for the class \code{VECM}.
#'
#' @param x An object of class \code{VECM}.
#'
#' @return A list with the following components: \cr
#'
#' \item{alpha}{The loading matrix.}
#' \item{beta}{The cointegrating matrix.}
#' \item{pi}{The long-run impact matrix.}
#' \item{gamma}{The short-run impact matrix.}
#' \item{D_ect}{The matrix containing the estimated deterministic terms within the error correction term.}
#' \item{D}{The matrix containing the estimated deterministic terms.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM}}
#'
#' @export
#'
#' @importFrom utils tail

coef.VECM <- function(x) {

  if (!inherits(x = x, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }

  K <- x$model$K
  p <- x$model$p

  beta <- x$beta
  beta_names <- rownames(beta)

  D_ect_coef <- NULL

  if (nrow(beta) > K) {
    beta_coef <- beta[1:K, ]
    if (is.vector(beta_coef)) {
      beta_coef <- as.matrix(beta_coef)
      colnames(beta_coef) <- beta_names[1]
    }
    D_ect_coef <- rbind(D_ect_coef, beta[-(1:K), ])
    rownames(D_ect_coef) <- beta_names[-(1:K)]
  } else {
    beta_coef <- beta
  }

  alpha <- x$alpha

  pi <- tcrossprod(x = alpha, y = beta_coef)

  gamma <- x$gamma
  gamma_names <- colnames(gamma)

  D_coef <- NULL

  if ("const" %in% gamma_names) {
    D_coef <- cbind(D_coef, "const" = gamma[, which("const" == gamma_names)])
    gamma <- gamma[, -which("const" == gamma_names)]
    gamma_names <- colnames(gamma)
  }

  if ("trend" %in% gamma_names) {
    D_coef <- cbind(D_coef, "trend" = gamma[, which("trend" == gamma_names)])
    gamma <- gamma[, -which("trend" == gamma_names)]
    gamma_names <- colnames(gamma)
  }

  season <- x$model$season

  if (!is.null(season)) {
    temp <- colnames(D_coef)
    temp_1 <- gamma_names[(K * p + 1):(K * p + season - 1)]
    D_coef <- cbind(D_coef, gamma[, (K * p + 1):(K * p + season - 1)])
    colnames(D_coef) <- c(temp, temp_1)
    gamma <- gamma[, -((K * p + 1):(K * p + season - 1))]
    gamma_names <- colnames(gamma)
  }

  exogen <- x$model$exogen

  if (!is.null(exogen)) {
    q <- x$model$q
    if (!is.null(q)) {
      temp <- ncol(gamma)
      D_coef <- cbind(D_coef, gamma[, (K * p + 1):temp])
      gamma <- gamma[, -((K * p + 1):temp)]
    } else {
      temp <- ncol(gamma)
      temp_1 <- colnames(D_coef)
      D_coef <- cbind(D_coef, gamma[, temp])
      colnames(D_coef) <- c(temp_1, tail(x = gamma_names, n = 1))
      gamma <- gamma[, -temp]
    }
  }

  return(list("alpha" = alpha, "beta" = beta_coef, "pi" = pi, "gamma" = gamma, "D_ect" = D_ect_coef, "D" = D_coef))

}
