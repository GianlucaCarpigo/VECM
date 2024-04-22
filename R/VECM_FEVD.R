#' Forecast error variance decomposition
#'
#' This function calculates the forecast error variance decomposition.
#'
#' @param object An object of class \code{VECM}.
#' @param horizon The time period over which the decomposition is performed.
#'
#' @return An object of class \code{VECM_FEVD} is a list with the following components: \cr
#'
#' \item{W}{The percentage of the forecast variance due to each random innovation.}
#' \item{forecast_error}{Forecast standard errors.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Lütkepohl, H., & Krätzig, M. (2004). \emph{Applied time series econometrics}. Cambridge: Cambridge University Press.
#'
#' Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
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
#' fevd <- VECM_FEVD(object = model, horizon = 10)
#' print(fevd)
#' }
#'
#' @seealso \code{\link{VECM}}
#'
#' @export

VECM_FEVD <- function(object, horizon = 10) {

  if (!inherits(x = object, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }

  VAR <- VECM_to_VAR(object)

  K <- object$model$K
  p <- object$model$p

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

  sigma_chol <- t(chol(sigma_u))

  theta <- array(data = NA, dim = c(K, K, horizon + 1))
  for (i in 1:(horizon + 1)) {
    theta[,, i] <- phi[,, i] %*% sigma_chol
  }

  B <- array(data = NA, dim = c(K, K, horizon + 1))
  B[,, 1] <- theta[,, 1]^2

  W_MSE <- matrix(data = NA, nrow = horizon + 1, ncol = K)
  W_MSE[1, ] <- rowSums(B[,, 1])

  W <- array(data = NA, dim = c(K, K, horizon + 1))
  W[,, 1] <- B[,, 1] / W_MSE[1, ]

  for (i in 2:(horizon + 1)) {
    B[,, i] <- B[,, i - 1] + theta[,, i]^2
    W_MSE[i, ] <- rowSums(B[,, i])
    W[,, i] <- B[,, i] / W_MSE[i, ]
  }

  var_names <- colnames(object$model$data)

  W <- aperm(a = W, perm = c(3, 2, 1))
  dimnames(W) <- list(NULL, var_names, var_names)

  colnames(W_MSE) <- var_names

  output <- list("W" = W[-(horizon + 1),,], "forecast_error" = sqrt(W_MSE[-(horizon + 1), ]))

  return(structure(.Data = output, class = "VECM_FEVD"))
}
