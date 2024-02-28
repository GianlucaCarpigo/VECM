#' Printing impulse response functions
#'
#' This function is a method for the class \code{VECM_IR}.
#'
#' @param x An object of class \code{VECM_IR}.
#' @param cum A logical value. If \code{TRUE} (not the default), accumulated responses are considered.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM_IR}}
#'
#' @export

print.VECM_IR <- function(x, cum = FALSE) {
  if (!inherits(x = x, what = "VECM_IR")) {
    stop("The object is not of class 'VECM_IR'.", call. = FALSE)
  }
  K <- dim(x$theta)[2]
  h <- dim(x$theta)[1]
  var_names <- colnames(x$theta[,, 1])
  method <- x$method
  transform <- x$transform
  if (method == "none") {
    cat("IMPULSE RESPONSE ANALYSIS\n\n")
    cat("transform:", transform, "\n\n")
    if (cum) {
      for(i in 1:K){
        cat("Response of", var_names[i], "\n\n")
        prmatrix(x = formatC(x = x$theta_cum[,, i], digits = 6, format = "f"), rowlab = paste0("h = ", 1:h), collab = var_names, quote = FALSE, right = TRUE)
        if (i != K) {
          cat("\n")
        }
      }
    } else {
      for(i in 1:K){
        cat("Response of", var_names[i], "\n\n")
        prmatrix(x = formatC(x = x$theta[,, i], digits = 6, format = "f"), rowlab = paste0("h = ", 1:h), collab = var_names, quote = FALSE, right = TRUE)
        if (i != K) {
          cat("\n")
        }
      }
    }
  } else {
    cat("IMPULSE RESPONSE ANALYSIS\n\n")
    cat("transform:", transform, "\n\n")
    cat("standard errors of the estimates are in parentheses\n\n")
    if (cum) {
      var_theta_cum <- matrix(data = 0, nrow = h, ncol = K^2)
      for (i in 2:h) {
        var_theta_cum[i, ] <- diag(x$sigma_theta_cum[,, i - 1])
      }
      temp <- array(data = var_theta_cum, dim = c(h, K, K))
      temp <- formatC(x = sqrt(temp), digits = 6, format = "f")
      mat <- matrix(data = NA, nrow = 2 * h, ncol = K)
      h_2 <- c()
      for (i in 1:K) {
        cat("Response of", var_names[i],"\n\n")
        for (j in 1:h) {
          mat[2 * j - 1, ] <- formatC(x = x$theta_cum[j,, i], digits = 6, format = "f")
          mat[2 * j, ] <- paste0("(", temp[j,, i], ")")
          h_2[2 * j - 1] <- paste0("h = ", j)
          h_2[2 * j] <- " "
        }
        prmatrix(x = mat, rowlab = h_2, collab = var_names, quote = FALSE, right = TRUE)
        if (i != K) {
          cat("\n")
        }
      }
    } else {
      if (transform == "none") {
        var_theta <- matrix(data = 0, nrow = h, ncol = K^2)
        for (i in 2:h) {
          var_theta[i, ] <- diag(x$sigma_theta[,, i - 1])
        }
      } else {
        var_theta <- matrix(data = NA, nrow = h, ncol = K^2)
        for (i in 1:h) {
          var_theta[i, ] <- diag(x$sigma_theta[,, i])
        }
      }
      temp <- array(data = var_theta, dim = c(h, K, K))
      temp <- formatC(x = sqrt(temp), digits = 6, format = "f")
      mat <- matrix(data = NA, nrow = 2 * h, ncol = K)
      h_2 <- c()
      for (i in 1:K) {
        cat("Response of", var_names[i],"\n\n")
        for (j in 1:h) {
          mat[2 * j - 1, ] <- formatC(x = x$theta[j,, i], digits = 6, format = "f")
          mat[2 * j, ] <- paste0("(", temp[j,, i], ")")
          h_2[2 * j - 1] <- paste0("h = ", j)
          h_2[2 * j] <- " "
        }
        prmatrix(x = mat, rowlab = h_2, collab = var_names, quote = FALSE, right = TRUE)
        if (i != K) {
          cat("\n")
        }
      }
    }
  }
}
