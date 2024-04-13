#' From VECM to VAR
#'
#' This function transforms a VEC model into a level VAR model.
#'
#' @param object An object of class \code{VECM}.
#'
#' @return An object of class \code{VAR} with the following components: \cr
#'
#' \item{A}{A matrix containing the VAR coefficients.}
#' \item{sigma_A}{The covariance matrix of the matrix \code{A}.}
#' \item{D}{A matrix containing the deterministic components.}
#' \item{Z}{A matrix containing both the current levels and the corresponding lagged levels data.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
#'
#' @seealso \code{\link{VECM}}
#'
#' @export

VECM_to_VAR <- function(object) {

  if (!inherits(x = object, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }

  K <- object$model$K
  p <- object$model$p

  model_coef <- coef.VECM(object)

  gamma <- array(data = model_coef$gamma, dim = c(K, K, p))

  pi <- model_coef$pi

  A <- array(data = NA, dim = c(K, K, p + 1))
  A[,, 1] <- pi + diag(x = 1, nrow = K, ncol = K) + gamma[,, 1]
  A[,, p + 1] <- -gamma[,, p]

  if (p > 1) {
    for (i in 2:p) {
      A[,, i] <- gamma[,, i] - gamma[,, i - 1]
    }
  }

  A <- matrix(data = A, nrow = K, ncol = K * (p + 1))

  data <- object$model$data

  var_names <- colnames(data)
  temp <- NULL
  for (i in 1:(p + 1)) {
    temp <- c(temp, paste0(var_names, "_L", i))
  }
  dimnames(A) <- list(var_names, temp)

  if (is.null(model_coef$D_ect)) {
    deter_ect <- NULL
  } else {
    deter_ect <- object$alpha %*% t(model_coef$D_ect)
    rownames(deter_ect) <- var_names
  }
 
  if (is.null(model_coef$D)) {
    deter <- NULL
  } else {
    deter <- model_coef$D
    rownames(deter) <- var_names
  }
  
  D <- cbind(deter, deter_ect)

  n <- object$model$n
  u <- object$u

  Z <- t(embed(x = data, dimension = p + 2)[, -(1:K)])

  temp1 <- NULL
  for (i in 1:K) {
    temp1 <- c(temp1, paste0(var_names[i], ":", temp))
  }

  sigma_A <- kronecker(X = solve(Z %*% t(Z)), Y = u %*% t(u))
  dimnames(sigma_A) <- list(temp1, temp1)
  output <- list("A" = A, "sigma_A" = sigma_A, "D" = D, "Z" = Z)

  return(structure(.Data = output, class = "VAR"))
}
