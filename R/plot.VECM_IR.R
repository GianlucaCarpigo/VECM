#' Plotting VEC model IRFs
#'
#' This function is a method for the class \code{VECM_IR}.
#' 
#' @param x An object of class \code{VECM_IR}.
#' @param type The type of plot. It can be chosen between \code{area} (the default), and \code{line}.
#' @param cum A logical value. If \code{TRUE} (not the default), accumulated responses are considered.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM_IR}}
#'
#' @export

plot.VECM_IR <- function(x, type = c("area", "line"), cum = FALSE) {
  if (!inherits(x = x, what = "VECM_IR")) {
    stop("The object is not of class 'VECM_IR'.", call. = FALSE)
  }
  K <- dim(x$theta)[2]
  h <- dim(x$theta)[1]
  x_seq <- 1:h
  var_names <- colnames(x$theta[,, 1])
  par(mfrow = c(K, K), mar = c(2, 2, 1.5, 1.5))
  if (x$method == "none") {
    if (cum) {
      for (i in 1:K) {
        for (j in 1:K) {
          temp <- paste0("Cumulative response of ", var_names[i], " to ", var_names[j])
          plot(x = x_seq, y = x$theta_cum[, j, i], xlab = "", ylab = "", main = temp, type = "l", lend = 1)
          title(xlab = "horizon", line = 2.4)
          abline(h = 0, lty = 2, lend = 1)
        }
      }
    } else {
      for (i in 1:K) {
        for (j in 1:K) {
          temp <- paste0("Response of ", var_names[i], " to ", var_names[j])
          plot(x = x_seq, y = x$theta[, j, i], xlab = "", ylab = "", main = temp, type = "l", lend = 1)
          title(xlab = "horizon", line = 2.4)
          abline(h = 0, lty = 2, lend = 1)
        }
      }
    }
  } else {
    if (cum) {
      var_theta_cum <- matrix(data = 0, nrow = h, ncol = K^2)
      for (i in 2:h) {
        var_theta_cum[i, ] <- diag(x$sigma_theta_cum[,, i - 1])
      }
      sd <- array(data = sqrt(var_theta_cum), dim = c(h, K, K))
    } else {
      if (x$transform == "none") {
        var_theta <- matrix(data = 0, nrow = h, ncol = K^2)
        for (i in 2:h) {
          var_theta[i, ] <- diag(x$sigma_theta[,, i - 1])
        }
      } else {
        var_theta <- matrix(data = NA, nrow = h, ncol = K^2)
        for (i in x_seq) {
          var_theta[i, ] <- diag(x$sigma_theta[,, i])
        }
      }
      sd <- array(data = sqrt(var_theta), dim = c(h, K, K))
    }
    alpha <- 0.05
    q_norm <- qnorm(c(alpha / 2, 1 - alpha / 2))
    lower <- array(data = NA, dim = dim(sd))
    dimnames(lower) <- list(NULL, var_names, var_names)
    upper <- array(data = NA, dim = dim(sd))
    dimnames(upper) <- list(NULL, var_names, var_names)
    type <- match.arg(type)
    if (cum) {
      for (i in 1:K) {
        for (j in 1:K) {
          lower[, j, i] <- x$theta_cum[, j, i] + q_norm[1] * sd[, j, i]
          upper[, j, i] <- x$theta_cum[, j, i] + q_norm[2] * sd[, j, i]
        }
      }
      if (type == "line") {
        for (i in 1:K) {
          for (j in 1:K) {
            temp <- paste0("Cumulative response of ", var_names[i], " to ", var_names[j])
            limits <- c(min(lower[, j, i]), max(upper[, j, i]))
            plot(x = x_seq, y = x$theta_cum[, j, i], ylim = limits, xlab = "", ylab = "", main = temp, type = "l", lend = 1)
            points(x = x_seq, y = lower[, j, i], type = "l", col = "red", lend = 1)
            points(x = x_seq, y = upper[, j, i], type = "l", col = "red", lend = 1)
            title(xlab = "horizon", line = 2.4)
            abline(h = 0, lty = 2, lend = 1)
          }
        }
      } else {
        x_lim <- c(x_seq, rev(x_seq))
        for (i in 1:K) {
          for (j in 1:K) {
            temp <- paste0("Cumulative response of ", var_names[i], " to ", var_names[j])
            limits <- c(min(lower[, j, i]), max(upper[, j, i]))
            plot(NA, xlim = c(1, h), ylim = limits, xlab = "", ylab = "", main = temp, cex.main = 0.9)
            title(xlab = "horizon", line = 2.4)
            y_1 <- lower[, j, i]
            y_2 <- upper[, j, i]
            polygon(x = x_lim, y = c(y_1, rev(y_2)), border = NA, col = "lightgrey")
            points(x = x_seq, y = x$theta_cum[, j, i], type = "l", lend = 1)
            abline(h = 0, lty = 2, lend = 1)
          }
        }
      }
    } else {
      for (i in 1:K) {
        for (j in 1:K) {
          lower[, j, i] <- x$theta[, j, i] + q_norm[1] * sd[, j, i]
          upper[, j, i] <- x$theta[, j, i] + q_norm[2] * sd[, j, i]
        }
      }
      if (type == "line") {
        for (i in 1:K) {
          for (j in 1:K) {
            temp <- paste0("Response of ", var_names[i], " to ", var_names[j])
            limits <- c(min(lower[, j, i]), max(upper[, j, i]))
            plot(x = x_seq, y = x$theta[, j, i], ylim = limits, xlab = "", ylab = "", main = temp, type = "l", lend = 1)
            points(x = x_seq, y = lower[, j, i], type = "l", col = "red", lend = 1)
            points(x = x_seq, y = upper[, j, i], type = "l", col = "red", lend = 1)
            title(xlab = "horizon", line = 2.4)
            abline(h = 0, lty = 2, lend = 1)
          }
        }
      } else {
        x_lim <- c(x_seq, rev(x_seq))
        for (i in 1:K) {
          for (j in 1:K) {
            temp <- paste0("Response of ", var_names[i], " to ", var_names[j])
            limits <- c(min(lower[, j, i]), max(upper[, j, i]))
            plot(NA, xlim = c(1, h), ylim = limits, xlab = "", ylab = "", main = temp, cex.main = 0.9)
            title(xlab = "horizon", line = 2.4)
            y_1 <- lower[, j, i]
            y_2 <- upper[, j, i]
            polygon(x = x_lim, y = c(y_1, rev(y_2)), border = NA, col = "lightgrey")
            points(x = x_seq, y = x$theta[, j, i], type = "l", lend = 1)
            abline(h = 0, lty = 2, lend = 1)
          }
        }
      }
    }
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}
