#' Printing VEC model forecasts
#'
#' This function is a method for the class \code{VECM_forecast}.
#'
#' @param x An object of class \code{VECM_forecast}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM_forecast}}
#'
#' @export

print.VECM_forecast <- function(x) {
  if (!inherits(x = x, what = "VECM_forecast")) {
    stop("The object is not of class 'VECM_forecast'.", call. = FALSE)
  }
  cat(paste0(toupper(x$type), " FORECAST VIA THE ESTIMATED VECTOR ERROR CORRECTION MODEL\n"))
  var_names <- colnames(x$data)
  u <- t(x$u)
  ME <- apply(X = u, MARGIN = 2, FUN = mean)
  RMSE <- sqrt(apply(X = u^2, MARGIN = 2, FUN = mean))
  MAE <- apply(X = abs(u), MARGIN = 2, FUN = mean)
  MPE <- apply(X = u / x$data, MARGIN = 2, FUN = mean) * 100
  MAPE <- apply(X = abs(u / x$data), MARGIN = 2, FUN = mean) * 100
  U_1 <- RMSE * 1 / (sqrt(apply(X = x$data^2, MARGIN = 2, FUN = mean)) + sqrt(apply(X = (x$data - u)^2, MARGIN = 2, FUN = mean)))
  U_2 <- sqrt(apply(X = (u[-1, ] / x$data[-nrow(x$data), ])^2, MARGIN = 2, FUN = mean) * 1 / apply(X = (diff(x$data) / x$data[-nrow(x$data), ])^2, MARGIN = 2, FUN = mean))
  if (x$type == "static") {
    sd <- diag(x$sigma_forecast)
    for (i in 1:length(var_names)) {
      cat(paste0("\nEquation ", i, ": ", var_names[i], "\n\n"))
      out <- cbind(x$forecast[, i], sd[i], x$CI_forecast[1, i], x$CI_forecast[2, i])
      prmatrix(x = formatC(x = out, digits = 5, format = "f"), rowlab = "T + 1 | T", collab = c("forecast", "se", "lower", "upper"), quote = FALSE, right = TRUE)
    }
    cat("\nFORECAST EVALUATION STATISTICS UP TO TIME T\n\n")
    out <- rbind(ME, RMSE, MAE, MPE, MAPE, U_1, U_2)
    prmatrix(x = formatC(x = out, digits = 5, format = "f"), rowlab = c("ME", "RMSE", "MAE", "MPE", "MAPE", "U1", "U2"), collab = var_names, quote = FALSE, right = TRUE)
  } else {
    row_names <- paste0("T + ", 1:nrow(x$forecast), " | T")
    sd <- aperm(a = sqrt(x$sigma_forecast), perm = c(3, 1, 2))
    for (i in 1:length(var_names)) {
      cat(paste0("\nEquation ", i, ": ", var_names[i], "\n\n"))
      out <- cbind(x$forecast[, i], sd[, i, i], x$CI_forecast[, 1, i], x$CI_forecast[, 2, i])
      prmatrix(x = formatC(x = out, digits = 5, format = "f"), rowlab = row_names, collab = c("forecast", "se", "lower", "upper"), quote = FALSE, right = TRUE)
    }
    cat("\nFORECAST EVALUATION STATISTICS UP TO TIME T\n\n")
    out <- rbind(ME, RMSE, MAE, MPE, MAPE, U_1, U_2)
    prmatrix(x = formatC(x = out, digits = 5, format = "f"), rowlab = c("ME", "RMSE", "MAE", "MPE", "MAPE", "U1", "U2"), collab = var_names, quote = FALSE, right = TRUE)
  }
}
