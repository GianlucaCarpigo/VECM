#' Printing cointegration test
#'
#' This function is a method for the class \code{cointegration_test}.
#'
#' @param x An object of class \code{cointegration_test}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{cointegration_test}}
#'
#' @export

print.cointegration_test <- function(x) {
  if (!inherits(x = x, what = "cointegration_test")) {
    stop("The object is not of class 'cointegration_test'.", call. = FALSE)
  }
  cat("COINTEGRATION TEST\n\n")
  method <- x$method
  if (method == "johansen") {
    cat("Method: Johansen\n")
  }
  if (method == "info_criteria") {
    cat("Method: information criteria\n")
  }
  cat("Model specification:", x$spec, "\n\n")
  cat("Number of observations:", x$n, "\n")
  cat("Lag order:", x$p, "\n\n")
  if (method == "johansen") {
    out <- formatC(x = cbind(x$e_val, x$stat, x$crit_val, x$p_val), digits = 5, format = "f")
    mat <- cbind("rank" = 0:(length(x$n_par) - 1), "param" = x$n_par, "eigenvalue" = out[, 1], "stat" = out[, 2], "quantile" = out[, 4], "p-value" = out[, 6])
    cat("TRACE TEST\n")
    prmatrix(x = mat, rowlab = rep(x = "", length = nrow(mat)), quote = FALSE, right = TRUE, na.print = " ")
    mat <- cbind("rank" = 0:(length(x$n_par) - 1), "param" = x$n_par, "eigenvalue" = out[, 1], "stat" = out[, 3], "quantile" = out[, 5], "p-value" = out[, 7])
    cat("\nMAXIMUM EIGENVALUE TEST\n")
    prmatrix(x = mat, rowlab = rep(x = "", length = nrow(mat)), quote = FALSE, right = TRUE, na.print = " ")
    cat("\n- quantiles are at the", 1 - x$alpha, "level\n")
    cat("- quantiles and p-values computed by sampling from the distributions of the test statistics\n")
    cat("- tests do not take exogenous variables into account")
  }
  if (method == "info_criteria") {
    out <- formatC(x = cbind(x$e_val, x$log_lik, x$IC), digits = 5, format = "f")
    mat <- cbind("rank" = 0:(length(x$n_par) - 1), "param" = x$n_par, "eigenvalue" = out[, 1], "log-lik" = out[, 2], out[, -(1:2)])
    prmatrix(x = mat, rowlab = rep(x = "", length = nrow(mat)), quote = FALSE, right = TRUE)
    pos <- apply(X = x$IC, MARGIN = 2, FUN = which.min) - 1
    cat(paste0("\n- chosen rank: r = ", pos[1], " (AIC) r = ", pos[2], " (BIC) r = ", pos[3], " (HQC)\n"))
    cat("- the test takes exogenous variables into account")
  }
}
