#' Plotting VEC model fits
#'
#' This function is a method for the class \code{VECM}.
#'
#' @param x An object of class \code{VECM}.
#' @param p_max The maximum lag used in the calculation of the (partial) autocorrelation function.
#' @param alpha The confidence level used in the calculation of the confidence intervals.
#' @param cross A logical value. If \code{TRUE} (not the default), cross correlations are considered.
#' @param kernel A character string giving the smoothing kernel used in the calculation of the Kernel density estimates. It can be selected between \code{biweight} (the default), \code{epanechnikov}, \code{gaussian}, \code{rectangular}, and \code{triangular}.
#'
#' @details The \code{plot.VECM} function displays different types of plots based on the user's selection. The options are the following:
#'
#'  1: error correction terms \cr
#'  2: raw residuals \cr
#'  3: raw residuals standardized \cr
#'  4: squared residuals \cr
#'  5: squared residuals standardized \cr
#'  6: autocorrelation functions of raw residuals \cr
#'  7: autocorrelation functions of squared residuals \cr
#'  8: partial autocorrelation functions of raw residuals \cr
#'  9: partial autocorrelation functions of squared residuals \cr
#' 10: kernel density estimations of raw residuals \cr
#'
#' In plotting the error correction term, you can further choose the following: \cr
#'
#' 1: representation without deterministic terms \cr
#' 2: representation with deterministic terms (if any) \cr
#' 3: representation without deterministic terms and without short-run dynamics \cr
#' 4: representation with deterministic terms (if any) but without short-run dynamics
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM}}
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
#' plot(model)
#'
#' }
#'
#' @export
#'
#' @importFrom graphics title
#' @importFrom stats acf density qnorm
#' @importFrom viridis viridis_pal

plot.VECM <- function(x, p_max = 25, alpha = 0.05, cross = FALSE, kernel = c("biweight", "epanechnikov", "gaussian", "rectangular", "triangular")) {

  if (!inherits(x = x, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }

  plot_type <- VECM_plot()
  plot_type <- match.arg(arg = plot_type, choices = c(as.character(1:10), ""))
  if (plot_type == "") {
    return()
  }

  if (plot_type == "1") {
    ECT_type <- ECT_plot()
    ECT_type <- match.arg(arg = ECT_type, choices = c(as.character(1:4), ""))
    if (ECT_type == "") {
      return()
    }
    K <- x$model$K
    N <- x$model$N
    p <- x$model$p
    r <- x$model$r
    if (ECT_type == "1") {
      ECT <- cbind(matrix(data = NA, nrow = r, ncol = p + 1), t(x$beta[1:K, ]) %*% x$model$LY[1:K, ])
      var_names <- colnames(x$beta)
      par(mfrow = c(1, 1))
      for (i in 1:r) {
        plot(x = 1:N, y = ECT[i, ], xlab = "", ylab = "", main = var_names[i], type = "l", lend = 1)
        title(xlab = "time", ylab = "error correction", line = 2.3)
        abline(h = 0, lty = 2, lend = 1)
      }
    }
    if (ECT_type == "2") {
      ECT <- cbind(matrix(data = NA, nrow = r, ncol = p + 1), t(x$beta) %*% x$model$LY)
      var_names <- colnames(x$beta)
      par(mfrow = c(1, 1))
      for (i in 1:r) {
        plot(x = 1:N, y = ECT[i, ], xlab = "", ylab = "", main = var_names[i], type = "l", lend = 1)
        title(xlab = "time", ylab = "error correction", line = 2.3)
        abline(h = 0, lty = 2, lend = 1)
      }
    }
    if (ECT_type == "3") {
      n <- x$model$n
      LDY <- x$model$LDY
      t_LDY <- t(LDY)
      M <- diag(x = 1, nrow = n, ncol = n) - t_LDY %*% solve(LDY %*% t_LDY) %*% LDY
      ECT <- cbind(matrix(data = NA, nrow = r, ncol = p + 1), t(x$beta[1:K, ]) %*% x$model$LY[1:K, ] %*% M)
      var_names <- colnames(x$beta)
      par(mfrow = c(1, 1))
      for (i in 1:r) {
        plot(x = 1:N, y = ECT[i, ], xlab = "", ylab = "", main = var_names[i], type = "l", lend = 1)
        title(xlab = "time", ylab = "error correction", line = 2.3)
        abline(h = 0, lty = 2, lend = 1)
      }
    }
    if (ECT_type == "4") {
      n <- x$model$n
      LDY <- x$model$LDY
      t_LDY <- t(LDY)
      M <- diag(x = 1, nrow = n, ncol = n) - t_LDY %*% solve(LDY %*% t_LDY) %*% LDY
      ECT <- cbind(matrix(data = NA, nrow = r, ncol = p + 1), t(x$beta) %*% x$model$LY %*% M)
      var_names <- colnames(x$beta)
      par(mfrow = c(1, 1))
      for (i in 1:r) {
        plot(x = 1:N, y = ECT[i, ], xlab = "", ylab = "", main = var_names[i], type = "l", lend = 1)
        title(xlab = "time", ylab = "error correction", line = 2.3)
        abline(h = 0, lty = 2, lend = 1)
      }
    }
  }

  if (plot_type == "2") {
    K <- x$model$K
    N <- x$model$N
    p <- x$model$p
    u <- cbind(matrix(data = NA, nrow = K, ncol = p + 1), x$u)
    var_names <- rownames(u)
    par(mfrow = c(1, 1))
    for (i in 1:K) {
      plot(x = 1:N, y = u[i, ], xlab = "", ylab = "", main = var_names[i], type = "l", lend = 1)
      title(xlab = "time", ylab = bquote(hat(u)), line = 2.3)
      abline(h = 0, lty = 2, lend = 1)
     }
  }

  if (plot_type == "3") {
    K <- x$model$K
    N <- x$model$N
    p <- x$model$p
    u <- x$u
    z <- cbind(matrix(data = NA, nrow = K, ncol = p + 1), (u - rowMeans(u)) / sqrt(diag(x$sigma_u)))
    q_norm <- qnorm(c(alpha / 2, 1 - alpha / 2))
    var_names <- rownames(u)
    par(mfrow = c(1, 1))
    for (i in 1:K) {
      y_lim <- c(min(q_norm[1], z[i, ], na.rm = TRUE), max(q_norm[2], z[i, ], na.rm = TRUE))
      plot(x = 1:N, y = z[i, ], ylim = y_lim, xlab = "", ylab = "", main = var_names[i], type = "l", lend = 1)
      title(xlab = "time", ylab = bquote(hat(z)), line = 2.3)
      abline(h = 0, lty = 2, lend = 1)
      abline(h = q_norm, lty = 2, col = "red", lend = 1)
    }
  }

  if (plot_type == "4") {
    K <- x$model$K
    N <- x$model$N
    p <- x$model$p
    u <- x$u
    u2 <- cbind(matrix(data = NA, nrow = K, ncol = p + 1), u^2)
    var_names <- rownames(u)
    par(mfrow = c(1, 1))
    for (i in 1:K) {
      plot(x = 1:N, y = u2[i, ], xlab = "", ylab = "", main = var_names[i], type = "l", lend = 1)
      title(xlab = "time", ylab = bquote(paste(hat(u), " squared")), line = 2.3)
      abline(h = 0, lty = 2, lend = 1)
    }
  }

  if (plot_type == "5") {
    K <- x$model$K
    N <- x$model$N
    p <- x$model$p
    u <- x$u
    z2 <- cbind(matrix(data = NA, nrow = K, ncol = p + 1), ((u - rowMeans(u)) / sqrt(diag(x$sigma_u)))^2)
    var_names <- rownames(u)
    par(mfrow = c(1, 1))
    for (i in 1:K) {
      plot(x = 1:N, y = z2[i, ], xlab = "", ylab = "", main = var_names[i], type = "l", lend = 1)
      title(xlab = "time", ylab = bquote(paste(hat(z), " squared")), line = 2.3)
      abline(h = 0, lty = 2, lend = 1)
    }
  }

  if (plot_type == "6") {
    K <- x$model$K
    n <- x$model$n
    u <- x$u
    ac_CI <- qnorm(c(alpha / 2, 1 - alpha / 2)) / sqrt(n)
    var_names <- rownames(u)
    if (cross) {
      ac <- acf(x = t(u), lag.max = p_max, type = "correlation", plot = FALSE)$acf
      mat_names <- matrix(data = NA, nrow = K, ncol = K)
      for (i in 1:K) {
        mat_names[i, ] <- paste0(var_names[i], " & ", var_names)
      }
      diag(mat_names) <- var_names
      for (i in 1:K) {
        par(mfrow = c(K, 1))
        for (j in 1:K) {
          y_lim = c(min(ac_CI[1], ac[, j, i]), max(ac_CI[2], ac[, j, i]))
          plot(x = 0:p_max, y = ac[, j, i], ylim = y_lim, xlab = "", ylab = "", main = mat_names[i, j], type = "h", lwd = 8, lend = 1)
          title(xlab = "lag", ylab =  bquote(paste("ac of ", hat(u))), line = 2.3)
          abline(h = 0, lty = 2, lend = 1)
          abline(h = ac_CI, lty = 2, col = "red", lend = 1)
        }
      }
    } else {
      par(mfrow = c(1, 1))
      for (i in 1:K) {
        ac <- acf(x = u[i, ], lag.max = p_max, type = "correlation", plot = FALSE)$acf
        y_lim = c(min(ac_CI[1], ac), max(ac_CI[2], ac))
        plot(x = 0:p_max, y = ac, ylim = y_lim, xlab = "", ylab = "", main = var_names[i], type = "h", lwd = 8, lend = 1)
        title(xlab = "lag", ylab =  bquote(paste("ac of ", hat(u))), line = 2.3)
        abline(h = 0, lty = 2, lend = 1)
        abline(h = ac_CI, lty = 2, col = "red", lend = 1)
      }
    }
  }

  if (plot_type == "7") {
    K <- x$model$K
    n <- x$model$n
    u <- x$u
    u2 <- u^2
    ac_CI <- qnorm(c(alpha / 2, 1 - alpha / 2)) / sqrt(n)
    var_names <- rownames(u)
    if (cross) {
      ac <- acf(x = t(u2), lag.max = p_max, type = "correlation", plot = FALSE)$acf
      mat_names <- matrix(data = NA, nrow = K, ncol = K)
      for (i in 1:K) {
        mat_names[i, ] <- paste0(var_names[i], " & ", var_names)
      }
      diag(mat_names) <- var_names
      for (i in 1:K) {
        par(mfrow = c(K, 1))
        for (j in 1:K) {
          y_lim = c(min(ac_CI[1], ac[, j, i]), max(ac_CI[2], ac[, j, i]))
          plot(x = 0:p_max, y = ac[, j, i], ylim = y_lim, xlab = "", ylab = "", main = mat_names[i, j], type = "h", lwd = 8, lend = 1)
          title(xlab = "lag", ylab =  bquote(paste("ac of ", hat(u), " squared")), line = 2.3)
          abline(h = 0, lty = 2, lend = 1)
          abline(h = ac_CI, lty = 2, col = "red", lend = 1)
        }
      }
    } else {
      par(mfrow = c(1, 1))
      for (i in 1:K) {
        ac <- acf(x = u2[i, ], lag.max = p_max, type = "correlation", plot = FALSE)$acf
        y_lim = c(min(ac_CI[1], ac), max(ac_CI[2], ac))
        plot(x = 0:p_max, y = ac, ylim = y_lim, xlab = "", ylab = "", main = var_names[i], type = "h", lwd = 8, lend = 1)
        title(xlab = "lag", ylab =  bquote(paste("ac of ", hat(u), " squared")), line = 2.3)
        abline(h = 0, lty = 2, lend = 1)
        abline(h = ac_CI, lty = 2, col = "red", lend = 1)
      }
    }
  }

  if (plot_type == "8") {
    K <- x$model$K
    n <- x$model$n
    u <- x$u
    pac_CI <- qnorm(c(alpha / 2, 1 - alpha / 2)) / sqrt(n)
    var_names <- rownames(u)
    if (cross) {
      pac <- acf(x = t(u), lag.max = p_max, type = "partial", plot = FALSE)$acf
      mat_names <- matrix(data = NA, nrow = K, ncol = K)
      for (i in 1:K) {
        mat_names[i, ] <- paste0(var_names[i], " & ", var_names)
      }
      diag(mat_names) <- var_names
      for (i in 1:K) {
        par(mfrow = c(K, 1))
        for (j in 1:K) {
          y_lim = c(min(pac_CI[1], pac[, j, i]), max(pac_CI[2], pac[, j, i]))
          plot(x = 0:p_max, y = c(0, pac[, j, i]), ylim = y_lim, xlab = "", ylab = "", main = mat_names[i, j], type = "h", lwd = 8, lend = 1)
          title(xlab = "lag", ylab =  bquote(paste("pac of ", hat(u))), line = 2.3)
          abline(h = 0, lty = 2, lend = 1)
          abline(h = pac_CI, lty = 2, col = "red", lend = 1)
        }
      }
    } else {
      par(mfrow = c(1, 1))
      for (i in 1:K) {
        pac <- acf(x = u[i, ], lag.max = p_max, type = "partial", plot = FALSE)
        y_lim = c(min(pac_CI[1], pac$acf), max(pac_CI[2], pac$acf))
        plot(x = 0:p_max, y = c(0, pac$acf), ylim = y_lim, xlab = "", ylab = "", main = var_names[i], type = "h", lwd = 8, lend = 1)
        title(xlab = "lag", ylab =  bquote(paste("pac of ", hat(u))), line = 2.3)
        abline(h = 0, lty = 2, lend = 1)
        abline(h = pac_CI, lty = 2, col = "red", lend = 1)
      }
    }
  }

  if (plot_type == "9") {
    K <- x$model$K
    n <- x$model$n
    u <- x$u
    u2 <- u^2
    pac_CI <- qnorm(c(alpha / 2, 1 - alpha / 2)) / sqrt(n)
    var_names <- rownames(u)
    if (cross) {
      pac <- acf(x = t(u2), lag.max = p_max, type = "partial", plot = FALSE)$acf
      mat_names <- matrix(data = NA, nrow = K, ncol = K)
      for (i in 1:K) {
        mat_names[i, ] <- paste0(var_names[i], " & ", var_names)
      }
      diag(mat_names) <- var_names
      for (i in 1:K) {
        par(mfrow = c(K, 1))
        for (j in 1:K) {
          y_lim = c(min(pac_CI[1], pac[, j, i]), max(pac_CI[2], pac[, j, i]))
          plot(x = 0:p_max, y = c(0, pac[, j, i]), ylim = y_lim, xlab = "", ylab = "", main = mat_names[i, j], type = "h", lwd = 8, lend = 1)
          title(xlab = "lag", ylab =  bquote(paste("pac of ", hat(u), " squared")), line = 2.3)
          abline(h = 0, lty = 2, lend = 1)
          abline(h = pac_CI, lty = 2, col = "red", lend = 1)
        }
      }
    } else {
      par(mfrow = c(1, 1))
      for (i in 1:K) {
        pac <- acf(x = u2[i, ], lag.max = p_max, type = "partial", plot = FALSE)
        y_lim = c(min(pac_CI[1], pac$acf), max(pac_CI[2], pac$acf))
        plot(x = 0:p_max, y = c(0, pac$acf), ylim = y_lim, xlab = "", ylab = "", main = var_names[i], type = "h", lwd = 8, lend = 1)
        title(xlab = "lag", ylab =  bquote(paste("pac of ", hat(u), " squared")), line = 2.3)
        abline(h = 0, lty = 2, lend = 1)
        abline(h = pac_CI, lty = 2, col = "red", lend = 1)
      }
    }
  }

  if (plot_type == "10") {
    kernel <- match.arg(kernel)
    K <- x$model$K
    u <- x$u
    var_names <- rownames(u)
    par(mfrow = c(1, 1))
    for (i in 1:K) {
      dens <- density(x = u[i, ], kernel = kernel)
      text <- paste0("h = ", round(x = dens$bw, digits = 4))
      plot(x = dens$x, y = dens$y, xlab = "", ylab = "", main = var_names[i], sub = text, type = "l", lend = 1)
      title(xlab = "values", line = 2.3)
      title(ylab = bquote(paste("density of  ", hat(u))), line = 2.3)
      abline(h = 0, lty = 2, lend = 1)
    }
  }
}
