#' Summarizing VEC model
#'
#' This function is a method for the class \code{VECM}.
#'
#' @param x An object of class \code{VECM}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM}}
#'
#' @export
#'
#' @importFrom stats printCoefmat symnum var
#' @importFrom matrixcalc upper.triangle

summary.VECM <- function(x) {

  if (!inherits(x = x, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }

  cat("-----------------------------\n")
  cat("VECTOR ERROR CORRECTION MODEL\n")
  cat("-----------------------------\n")

  method <- x$model$method

  cat("\nEstimation method:", method)

  spec <- x$model$spec

  cat("\nDeterministic terms:", spec)

  n <- x$model$n

  cat("\n\nNumber of observations after adjustments:", n)

  n_par <- x$model$n_par

  cat("\nNumber of estimated parameters:", n_par)

  p <- x$model$p

  cat("\n\nLag order:", p)

  r <- x$model$r

  cat("\nCointegration rank:", r)

  K <- x$model$K

  u <- x$u

  sigma_u <- x$sigma_u

  det_sigma_u <- det(sigma_u)
  log_lik <- -0.5 * (n * K * log(2 * pi) + n * log(det_sigma_u) + sum(diag(t(u) %*% solve(sigma_u) %*% u)))
  cat("\n\nlog-likelihood:", formatC(x = log_lik, digits = 5, format = "f"))
  cat("\ndet(sigma_u):", formatC(x = det_sigma_u, digits = 5, format = "e"))
  model_fit <- -2 * log_lik / n

  model_par <- n_par + r^2
  cat("\n\nAIC:", formatC(x = model_fit + 2 * model_par / n, digits = 5, format = "f"))
  cat("\nBIC:", formatC(x = model_fit + log(n) * model_par / n, digits = 5, format = "f"))
  cat("\nHQC:", formatC(x = model_fit + 2 * log(log(n)) * model_par / n, digits = 5, format = "f"))

  cat("\n\nSignificance levels: ***", 0.01, "**", 0.05, "*", 0.1)

  cat("\n\nCOINTEGRATING VECTORS\n")

  beta <- x$beta
  beta_names <- rownames(beta)

  beta_k_r <- beta[-(1:r), ]

  if (is.vector(beta_k_r)) {
    beta_k_r <- as.matrix(beta_k_r)
    rownames(beta_k_r) <- beta_names[-(1:r)]
  }

  sigma_beta <- x$sigma_beta

  beta_se <- matrix(data = sqrt(diag(sigma_beta)), nrow = nrow(beta_k_r), ncol = r, byrow = TRUE)

  beta_tval <- beta_k_r / beta_se

  d <- x$model$d

  beta_pval <- 2 * pt(q = abs(beta_tval), df = n - d, lower.tail = FALSE)

  coef_star <- apply(X = beta_pval, MARGIN = 2, FUN = symnum, cutpoints = c(0, 0.01, 0.05, 0.1, 1), symbols = c("***", "** ", "*  ", ""), legend = FALSE, numeric.x = TRUE, corr = FALSE)

  for (i in 1:r) {
    cat("\nNormalization with respect to:", beta_names[i], "\n\n")
    temp <- cbind(beta_k_r[, i], beta_se[, i], beta_tval[, i], beta_pval[, i])
    temp <- formatC(x = temp, digits = 5, format = "f")
    prmatrix(x = cbind(temp, coef_star[, i]), collab = c("estimate", "s.e.", "t-value", "p-value", ""), quote = FALSE, right = TRUE)
  }

  cat("\nDIFFERENCED VAR EQUATIONS\n")

  alpha <- x$alpha
  gamma <- x$gamma

  alphagamma <- t(cbind(alpha, gamma))
  alphagamma_row <- nrow(alphagamma)

  sigma_alphagamma <- x$sigma_alphagamma

  alphagamma_se <- matrix(data = sqrt(diag(sigma_alphagamma)), nrow = alphagamma_row, ncol = K, byrow = TRUE)

  alphagamma_tval <- alphagamma / alphagamma_se

  alphagamma_pval <- 2 * pt(q = abs(alphagamma_tval), df = n - d, lower.tail = FALSE)

  coef_star <- apply(X = alphagamma_pval, MARGIN = 2, FUN = symnum, cutpoints = c(0, 0.01, 0.05, 0.1, 1), symbols = c("***", "** ", "*  ", ""), legend = FALSE, numeric.x = TRUE, corr = FALSE)

  DY <- x$model$DY

  DY_mean <- rowMeans(DY)
  DY_no_mean <- DY - DY_mean

  DY_mean <- formatC(x = DY_mean, digits = 5, format = "f")

  DY_sd <- sqrt(apply(X = DY, MARGIN = 1, FUN = var))
  DY_sd <- formatC(x = DY_sd, digits = 5, format = "f")

  TSS <- diag(DY_no_mean %*% t(DY_no_mean))
  RSS <- diag(u %*% t(u))

  RSE <- sqrt(RSS / (n - alphagamma_row))
  RSE <- formatC(x = RSE, digits = 5, format = "f")

  R2 <- 1 - RSS / TSS
  R2_adj <- 1 - (1 - R2) * (n - 1) / (n - alphagamma_row)

  R2 <- formatC(x = R2, digits = 5, format = "f")
  R2_adj <- formatC(x = R2_adj, digits = 5, format = "f")

  DW_stat <- colSums(apply(X = u, MARGIN = 1, FUN = diff)^2) / RSS
  DW_stat <- formatC(x = DW_stat, digits = 5, format = "f")

  alpha_names <- rownames(alpha)

  for (i in 1:K) {
    cat("\nEquation ", i, ": ", alpha_names[i], "\n\n", sep = "")
    temp <- cbind(alphagamma[, i], alphagamma_se[, i], alphagamma_tval[, i], alphagamma_pval[, i])
    temp <- formatC(x = temp, digits = 5, format = "f")
    prmatrix(x = cbind(temp, coef_star[, i]), collab = c("estimate", "s.e.", "t-value", "p-value", ""), quote = FALSE, right = TRUE)
    cat("---")
    cat("\nE(", alpha_names[i], "): ", DY_mean[i], sep = "")
    cat("\nSD(", alpha_names[i], "): ", DY_sd[i], sep = "")
    cat("\nMultiple R-squared:", R2[i])
    cat("\nAdjusted R-squared:", R2_adj[i])
    cat("\nResidual standard error:", RSE[i])
    cat("\nDurbin-Watson statistic:", DW_stat[i])
    cat("\n")
  }

  cat("\nF-TEST TABLE (H0: all parameters = 0)\n\n")

  F_stat <- ((TSS - RSS) / (alphagamma_row - 1)) / (RSS / (n - alphagamma_row))
  F_pval <- pf(q = F_stat, df1 = alphagamma_row - 1, df2 = n - alphagamma_row, lower.tail = FALSE)

  coef_star <- symnum(x = F_pval, cutpoints = c(0, 0.01, 0.05, 0.1, 1), symbols = c("***", "** ", "*  ", ""), legend = FALSE, numeric.x = TRUE, corr = FALSE)

  temp <- cbind(TSS, RSS, F_stat, F_pval)
  temp <- formatC(x = temp, digits = 5, format = "f")

  prmatrix(x = cbind(temp, coef_star), rowlab = paste0("Equation ", 1:K), collab = c("TSS", "RSS", "F-value", "p-value", ""), quote = FALSE, right = TRUE)

  cat("---")
  cat("\nDF numerator:", alphagamma_row - 1)
  cat("\nDF denominator:", n - alphagamma_row)

  cat("\n\nRESIDUALS OF DIFFERENCED VAR EQUATIONS\n\n")

  temp <- t(apply(X = u, MARGIN = 1, FUN = summary))
  temp <- formatC(x = temp, digits = 5, format = "f")

  prmatrix(x = temp, rowlab = paste0("Equation ", 1:K), collab = c("minimum", "Q1", "median", "mean", "Q3", "maximum"), quote = FALSE, right = TRUE)

  cat("\nRESIDUAL COVARIANCE MATRIX BETWEEN DIFFERENCED VAR EQUATIONS\n\n")

  sigma_u <- formatC(x = sigma_u, digits = 5, format = "f")
  sigma_u <- upper.triangle(sigma_u)
  sigma_u[sigma_u == "0"] <- NA

  prmatrix(x = sigma_u, quote = FALSE, right = TRUE, na.print = "")

  cat("\nMOVING-AVERAGE IMPACT MATRIX\n\n")

  C <- VECM_MA(x)$C
  C <- formatC(x = C, digits = 5, format = "f")

  prmatrix(x = C, quote = FALSE, right = TRUE)

}
