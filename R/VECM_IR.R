#' Impulse responses analysis
#' 
#' This function calculates impulse response functions for VEC models.
#'
#' @param object An object of class \code{VECM}.
#' @param horizon The forecasting horizon.
#' @param transform ...
#' @param method ...
#'
#' @return An object of class \code{VECM_IR} is a list with the following components: \cr
#' 
#' \item{theta}{...} 
#' \item{sigma_theta}{...} 
#' \item{theta_cum}{...} 
#' \item{sigma_theta_cum}{...} 
#' \item{transform}{...} 
#' \item{method}{...}
#' 
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Lütkepohl, H., & Krätzig, M. (2004). \emph{Applied time series econometrics}. Cambridge: Cambridge University Press.
#'
#' Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
#' 
#' @export
#'
#' @examples #ciao
#'
#' @importFrom expm "%^%"
#' @importFrom matrixcalc D.matrix K.matrix L.matrix

VECM_IR <- function(object, horizon = 10, transform = c("none", "orthogonal"), method = c("none", "analytic")) {
  
  if (!inherits(x = object, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }
  
  VAR <- VECM_to_VAR(object)
  
  K <- object$model$K
  p <- object$model$p
  n <- object$model$n
  
  sigma_u <- object$sigma_u
  
  A <- array(data = VAR$A, dim = c(K, K, p + 1))
  
  A_new <- array(data = 0, dim = c(K, K, horizon))
  if (horizon < (p + 1)) {
    A_new <- A[,, 1:horizon]
  } else {
    A_new[,, 1:(p + 1)] <- A
  }
  
  id_mat <- diag(x = 1, nrow = K, ncol = K)
  
  phi <- array(data = 0, dim = c(K, K, horizon + 1))
  phi[,, 1] <- id_mat
  for (i in 2:(horizon + 1)) {
    for (j in 1:(i - 1)) {
      phi[,, i] <- phi[,, i] + phi[,, i - j] %*% A_new[,, j]
    }
  }
  
  phi_cum <- array(data = NA, dim = c(K, K, horizon + 1))
  phi_cum[,, 1] <- id_mat
  for (i in 2:(horizon + 1)) {
    phi_cum[,, i] <- phi_cum[,, i - 1] + phi[,, i]
  }
  
  transform <- match.arg(transform)
  
  if (transform == "none") {
    theta <- phi
    theta_cum <- phi_cum
  }
  
  if (transform == "orthogonal") {
    sigma_chol <- t(chol(sigma_u))
    theta <- array(data = NA, dim = c(K, K, horizon + 1))
    theta_cum <- array(data = NA, dim = c(K, K, horizon + 1))
    for (i in 1:(horizon + 1)) {
      theta[,, i] <- phi[,, i] %*% sigma_chol
      theta_cum[,, i] <- phi_cum[,, i] %*% sigma_chol
    }
  }
  
  method <- match.arg(method)
  
  if (method == "none") {
    sigma_theta <- NULL
    sigma_theta_cum <- NULL
  }
  
  if (method == "analytic") {
    A <- VAR$A
    sigma_A <- VAR$sigma_A
    R <- rbind(A, cbind(diag(x = 1, nrow = K * p, ncol = K * p), matrix(data = 0, nrow = K * p, ncol = K)))
    t_R <- t(R)
    J <- cbind(id_mat, matrix(data = 0, nrow = K, ncol = K * p))
    G <- array(data = NA, dim = c(K^2, (p + 1) * K^2, horizon))
    for (i in 1:horizon) {
      temp <- matrix(data = 0, nrow = K^2, ncol = (p + 1) * K^2)
      index <- 0:(i - 1)
      for (j in 1:length(index)) {
        temp <- temp + kronecker(X = J %*% (t_R %^% (i - 1 - index[j])), Y = phi[,, j])
      }
      G[,, i] <- temp
    }
    G_cum <- array(data = NA, dim = c(K^2, (p + 1) * K^2, horizon))
    G_cum[,, 1] <- G[,, 1]
    for (i in 2:horizon) {
      G_cum[,, i] <- G_cum[,, i - 1] + G[,, i]
    }
    
    if (transform == "none") {
      sigma_theta <- array(data = NA, dim = c(K^2, K^2, horizon))
      sigma_theta_cum <- array(data = NA, dim = c(K^2, K^2, horizon))
      for (i in 1:horizon) {
        sigma_theta[,, i] <- G[,, i] %*% sigma_A %*% t(G[,, i])
        sigma_theta_cum[,, i] <- G_cum[,, i] %*% sigma_A %*% t(G_cum[,, i])
      }
    }
    
    if (transform == "orthogonal") {
      C <- array(data = NA, dim = c(K^2, (p + 1) * K^2, horizon + 1))
      C[,, 1] <- 0
      C_bar <- array(data = NA, dim = c(K^2, K * (K + 1) / 2, horizon + 1))
      L_mat <- L.matrix(K)
      t_L_mat <- t(L_mat)
      H <- t_L_mat %*% solve(L_mat %*% (diag(x = 1, nrow = K^2, ncol = K^2) + K.matrix(r = K, c = K)) %*% kronecker(X = sigma_chol, Y = id_mat) %*% t_L_mat)
      C_bar[,, 1] <- kronecker(X = id_mat, Y = phi[,, 1]) %*% H
      t_sigma_chol <- t(sigma_chol)
      for (i in 2:(horizon + 1)) {
        C[,, i] <- kronecker(X = t_sigma_chol, Y = id_mat) %*% G[,, i - 1]
        C_bar[,, i] <- kronecker(X = id_mat, Y = phi[,, i]) %*% H
      }
      D_mat <- D.matrix(K)
      t_D_mat <- t(D_mat)
      D_mat_plus <- solve(t_D_mat %*% D_mat) %*% t_D_mat
      sigma_s <- 2 * D_mat_plus %*% kronecker(X = sigma_u, Y = sigma_u) %*% t(D_mat_plus)
      B <- array(data = NA, dim = c(K^2, (p + 1) * K^2, horizon))
      B_bar <- array(data = NA, dim = c(K^2, K * (K + 1) / 2, horizon))
      for (i in 1:horizon) {
        B[,, i] <- kronecker(X = t_sigma_chol, Y = id_mat) %*% G_cum[,, i]
        B_bar[,, i] <- kronecker(X = id_mat, Y = phi_cum[,, i + 1]) %*% H
      }
      sigma_theta <- array(data = NA, dim = c(K^2, K^2, horizon + 1))
      for (i in 1:(horizon + 1)) {
        sigma_theta[,, i] <- C[,, i] %*% sigma_A %*% t(C[,, i]) + C_bar[,, i] %*% sigma_s %*% t(C_bar[,, i])
      }
      sigma_theta_cum <- array(data = NA, dim = c(K^2, K^2, horizon))
      for (i in 1:horizon) {
        sigma_theta_cum[,, i] <- B[,, i] %*% sigma_A %*% t(B[,, i]) + B_bar[,, i] %*% sigma_s %*% t(B_bar[,, i])
      }
    }
  }
  
  theta <- aperm(a = theta, perm = c(3, 1, 2))
  theta_cum <- aperm(a = theta_cum, perm = c(3, 1, 2))
  
  var_names <- colnames(object$model$data)
  dimnames(theta) <- list(NULL, var_names, var_names)
  dimnames(theta_cum) <- list(NULL, var_names, var_names)
  
  output <- list("theta" = theta, "sigma_theta" = sigma_theta / n, "theta_cum" = theta_cum, "sigma_theta_cum" = sigma_theta_cum / n, "transform" = transform, "method" = method)
  
  return(structure(.Data = output, class = "VECM_IR"))
  
}
