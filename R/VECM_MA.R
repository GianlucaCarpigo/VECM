#' Moving average impact matrix
#'
#' This function provides the MA representation of the VEC model (i.e. Granger-representation theorem).
#'
#' @param object An object of class \code{VECM}.
#'
#' @return A list containing the following components: \cr
#'
#' \item{C}{The moving average impact matrix.}
#' \item{sigma_C}{The covariance matrix of the matrix \code{C}.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM}}
#'
#' @export

VECM_MA <- function(object) {

  if (!inherits(x = object, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }

  K <- object$model$K

  id_mat <- diag(x = 1, nrow = K, ncol = K)

  model_coef <- coef(object)
  
  gamma <- model_coef$gamma

  p <- object$model$p

  index <- seq(from = 1, to = K * p, by = K)

  temp <- matrix(data = 0, nrow = K, ncol = K)

  for (i in 1:length(index)) {
    temp <- temp + gamma[, index[i]:(index[i] + K - 1)]
  }

  gamma_hat <- id_mat - temp

  r <- object$model$r

  alpha <- model_coef$alpha
  alpha_oc <- orth_compl(mat = alpha, r = r)

  beta <- model_coef$beta
  beta_oc <- orth_compl(mat = beta, r = r)

  C <- beta_oc %*% solve(t(alpha_oc) %*% gamma_hat %*% beta_oc) %*% t(alpha_oc)

  xi_1 <- (t(gamma_hat %*% C) - id_mat) %*% alpha %*% solve(t(alpha) %*% alpha)
  xi_2 <- kronecker(X = matrix(data = 1, nrow = p, ncol = 1), Y = C, make.dimnames = TRUE)
  xi <- cbind(xi_1, t(xi_2))

  S_q <- kronecker(X = C, Y = xi, make.dimnames = TRUE)

  sigma_alphagamma <- object$sigma_alphagamma

  spec <- object$model$spec

  if (spec == "uconst" | spec == "rtrend") {
    sigma_names <- rownames(sigma_alphagamma)
    temp <- paste0("const:", rownames(alpha))
    temp <- which(temp %in% sigma_names)
    sigma_alphagamma <- sigma_alphagamma[-temp, -temp]
  }

  if (spec == "utrend") {
    sigma_names <- rownames(sigma_alphagamma)
    alpha_names <- rownames(alpha)
    temp <- c(paste0("const:", alpha_names), paste0("trend:", alpha_names))
    temp <- which(temp %in% sigma_names)
    sigma_alphagamma <- sigma_alphagamma[-temp, -temp]
  } 
  
  season <- object$model$season
  
  if (!is.null(season)) {
    sigma_names <- rownames(sigma_alphagamma)
    alpha_names <- rownames(alpha)
    temp <- NULL
    for (i in 1:K) {
      temp <- c(temp, paste0("season_", c(1:season), ":", alpha_names[i]))
    }
    temp <- which(temp %in% sigma_names)
    sigma_alphagamma <- sigma_alphagamma[-temp, -temp]
  }
  
  exogen_names <- rownames(object$model$exogen)
  
  if (!is.null(exogen_names)){
    sigma_names <- rownames(sigma_alphagamma)
    alpha_names <- rownames(alpha)
    temp <- NULL
    for (i in 1:K) {
      temp <- c(temp, paste0(exogen_names, ":", alpha_names[i]))
    }
    temp <- which(temp %in% sigma_names)
    sigma_alphagamma <- sigma_alphagamma[-temp, -temp]
  }

  d <- object$model$d
  n <- object$model$n

  sigma_C <- (n - d) / n * S_q %*% sigma_alphagamma %*% t(S_q)

  return(list("C" = C, "sigma_C" = sigma_C))

}
