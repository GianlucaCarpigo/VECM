#' Forecasting using VEC models
#'
#' This function returns forecasts for VEC models.
#'
#' @param object An object of class \code{VECM}.
#' @param horizon The forecasting horizon.
#' @param type The type of forecast. It can be chosen between \code{static} (the default) and \code{dynamic}.
#' @param alpha The confidence level of the forecast confidence intervals.
#'
#' @return An object of class \code{VECM_forecast} is a list with the following components: \cr
#'
#' \item{forecast}{The matrix containing the forecast of the levels variables.}
#' \item{sigma_forecast}{The covariance matrix of the matrix \code{forecast}.}
#' \item{CI_forecast}{The array containing the confidence intervals of the forecasts.}
#' \item{type}{The type of forecast.}
#' \item{alpha}{The confidence level of the forecast confidence intervals.}
#' \item{data}{The matrix containing the observed values.}
#' \item{u}{The residuals matrix.}
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
#' #' require(bvartools)
#'
#' ## The data for this example are taken from Lütkepohl, 2005.
#' data("e6")
#'
#' model <- VECM(data = e6, p = 3,  r = 1, method = "EGLS", spec = "uconst", season = 4)
#'
#' forecast <- VECM_forecast(object = model, type = "dynamic")
#'
#' print(forecast)
#'
#' plot(forecast)
#' }
#'
#' @export
#'
#' @importFrom expm "%^%"
#' @importFrom stats qnorm

VECM_forecast <- function(object, horizon = 10, type = c("static", "dynamic"), alpha = 0.05) {

  if (!inherits(x = object, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }

  VAR <- VECM_to_VAR(object)

  A <- VAR$A

  K <- object$model$K
  p <- object$model$p
  n <- object$model$n

  data <- object$model$data[-(1:(p + 1)), ]

  sigma_u <- object$sigma_u

  type <- match.arg(type)

  if (type == "static") {
    if (horizon > 1) {
      warning(paste0("For 'type' = '", type, "' forecast is only one-step ahead."), call. = FALSE)
    }
    y <- apply(X = data, MARGIN = 2, FUN = rev)[1:(p + 1), ]
    y <- c(t(y))
    data_new <- t(A %*% y)
    sigma_data_new <- (n + K * p + 1) / n * sigma_u
    sd <- sqrt(diag(sigma_data_new))
    CI <- matrix(data = NA, nrow = 2, ncol = K)
    dimnames(CI) <- c(list(c("lower", "upper"), colnames(data)))
    q_norm <- qnorm(c(alpha / 2, 1 - alpha / 2))
    for (i in 1:K) {
      CI[, i] <- c(data_new[, i] + q_norm[1] * sd[i], data_new[, i] + q_norm[2] * sd[i])
    }
  }

  if (type == "dynamic") {
    if (horizon == 1) {
      warning("Since 'horizon' = 1 you are not doing a dynamic forecast. Please select 'type' = 'static'.", call. = FALSE)
    }
    data_new <- matrix(data = NA, nrow = horizon, ncol = K)
    colnames(data_new) <- colnames(data)
    for (i in 1:horizon) {
      y <- apply(X = data, MARGIN = 2, FUN = rev)[1:(p + 1), ]
      y <- c(t(y))
      data_new[i, ] <- t(A %*% y)
      data <- rbind(data, data_new[i, ])
    }
    B <- rbind(A, cbind(diag(x = 1, nrow = K * p, ncol = K * p), matrix(data = 0, nrow = K * p, ncol = K)))
    A <- array(data = A, dim = c(K, K, p + 1))
    A_new <- array(data = 0, dim = c(K, K, horizon))
    if (horizon < (p + 1)) {
      A_new <- A[,, 1:horizon]
    } else {
      A_new[,, 1:(p + 1)] <- A
    }
    phi <- array(data = 0, dim = c(K, K, horizon + 1))
    phi[,, 1] <- diag(x = 1, nrow = K, ncol = K)
    for (i in 2:(horizon + 1)) {
      for (j in 1:(i - 1)) {
        phi[,, i] <- phi[,, i] + phi[,, i - j] %*% A_new[,, j]
      }
    }
    MSE <- array(data = NA, dim = c(K, K, horizon))
    MSE[,, 1] <- phi[,, 1] %*% sigma_u %*% t(phi[,, 1])
    for (i in 2:horizon) {
      MSE[,, i] <- MSE[,, i - 1] + phi[,, i] %*% sigma_u %*% t(phi[,, i])
    }
    J <- cbind(0, diag(x = 1, nrow = K, ncol = K), matrix(data = 0, nrow = K, ncol = K * p))[, -(K * (p + 1) + 1)]
    t_J <- t(J)
    t_B <- t(B)
    Z <- VAR$Z
    gamma <- Z %*% t(Z)
    gamma_inv <- solve(gamma)
    omega <- array(data = NA, dim = c(K, K, horizon))
    for (h in 1:horizon) {
      temp <- matrix(data = 0, nrow = K, ncol = K)
      index <- 0:(h - 1)
      for (i in 1:length(index)) {
        for (j in 1:length(index)) {
          temp <- temp + sum(diag((t_B %^% (h - 1 - index[i])) %*% gamma_inv %*% (B %^% (h - 1 - index[j])) %*% gamma)) * (J %*% (B %^% index[i]) %*% t_J) %*% sigma_u %*% t(J %*% (B %^% index[j]) %*% t_J)
        }
      }
      omega[,, h] <- temp
    }
    sigma_data_new <- array(data = NA, dim = c(K, K, horizon))
    dimnames(sigma_data_new) <- list(colnames(data), colnames(data), NULL)
    sd <- matrix(data = NA, nrow = horizon, ncol = K)
    for (i in 1:horizon) {
      sigma_data_new[,, i] <- MSE[,, i] + omega[,, i] / n
      sd[i, ] <- sqrt(diag(sigma_data_new[,, i]))
    }
    CI <- array(data = NA, dim = c(horizon, 2, K))
    dimnames(CI) <- list(NULL, c("lower", "upper"), colnames(data))
    q_norm <- qnorm(c(alpha / 2, 1 - alpha / 2))
    for (i in 1:K) {
      for (j in 1:horizon) {
        CI[j, , i] <- c(data_new[j, i] + q_norm[1] * sd[j, i], data_new[j, i] + q_norm[2] * sd[j, i])
      }
    }
  }

  output <- list("forecast" = data_new, "sigma_forecast" = sigma_data_new, "CI_forecast" = CI, "type" = type, "alpha" = alpha, "data" = object$model$data[-(1:(p + 1)), ], "u" = object$u)

  return(structure(.Data = output, class = "VECM_forecast"))

}
