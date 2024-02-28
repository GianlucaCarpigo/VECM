#' Cointegration test
#'
#' This function implements a test to estimate the cointegrating rank between variables. The user can choose between the Johansen test and an information criteria-based test.
#'
#' @param data A matrix containing the multivariate time series to be analyzed.
#' @param p The lag order.
#' @param spec The type of deterministic component to be considered. It can be selected between \code{none} (the default), \code{rconst}, \code{uconst}, \code{rtrend}, and \code{utrend}.
#' @param method The estimation method of the cointegrating rank. It can be selected between \code{johansen} (the default) and \code{info_criteria}.
#' @param alpha The significance level used for hypothesis testing.
#' @param shift An integer or a vector. The break point at which there is a change in the levels of the variables.
#' @param season The number of seasonal dummy variables.
#' @param season_ref The season selected as reference. If \code{NULL}, the last season is automatically used as reference.
#' @param exogen_ect The exogenous variable to be included in the error correction term.
#' @param exogen The exogenous variable to be included outside the error correction term.
#'
#' @details When choosing \code{method} = "johansen", the function incorporates both the trace and maximum eigenvalues tests. The determination of the cointegrating rank should consider the results from both tests.
#'
#' When choosing \code{method} = "info_criteria", the function computes the most used information criteria, i.e. AIC, BIC, and HQC.
#'
#' @return An object of class \code{cointegration_test} is a list with the following components: \cr
#'
#' \item{n}{The number of observations.}
#' \item{p}{The lag order.}
#' \item{spec}{The type of deterministic component.}
#' \item{method}{The estimation method of the cointegrating rank.}
#' \item{n_par}{A vector containing the number of estimated parameters.}
#' \item{e_val}{A vector cotaining the eigenvalues used to implement the test.}
#' \item{stat}{A matrix containing the Johansen test statitics.}
#' \item{crit_val}{A matrix containing the critical values associated with the Johansen distribution.}
#' \item{p_val}{A matrix containing the p-values of the Johansen test statistics.}
#' \item{alpha}{The probability level at which the critical values are computed.}
#' \item{log_lik}{A vector containing the log-likelihood of the estimated VEC models.}
#' \item{IC}{A matrix containing the information criteria based on \code{log_lik}.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Johansen, S. (1995). \emph{Likelihood-based Inference in Cointegrated Vector Autoregressive Models}. Oxford: Oxford University Press.
#'
#' Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
#'
#' @seealso \code{\link{johansen_distr}} \code{\link{VECM}}
#'
#' @examples
#' \dontrun{
#' require(bvartools)
#'
#' ## The data for this example are taken from Lütkepohl, 2005.
#' data("e6")
#'
#' co_test <- cointegration_test(data = e6, p = 3, spec = "uconst", method = "johansen")
#'
#' print(co_test)
#'
#' }
#'
#' @export
#'
#' @importFrom stats quantile

cointegration_test <- function(data, p, spec = c("none", "rconst", "uconst", "rtrend", "utrend"), method = c("johansen", "info_criteria"), alpha = 0.05, shift = NULL, season = NULL, season_ref = NULL, exogen_ect = NULL, exogen = NULL) {

  data <- as.matrix(data)

  if (NA %in% data) {
    stop("'data' contains missing values.", call. = FALSE)
  }

  K <- ncol(data)
  N <- nrow(data)
  n <- N - (p + 1)

  LY <- t(data[(1 + p):(N - 1), ])
  Z <- embed(x = diff(data), dimension = p + 1)
  DY <- t(Z[, 1:K])
  LDY <- t(Z[, -(1:K)])

  r <- 0:K

  spec <- match.arg(spec)

  if (spec == "none") {
    n_par <- K * r + K * r + K^2 * p
  }

  if (spec == "rconst") {
    LY <- rbind(LY, "const" = 1)
    n_par <- K * r + (K + 1) * r + K^2 * p
  }

  if (spec == "uconst") {
    LDY <- rbind(LDY, "const" = 1)
    n_par <- K * r + K * r + K * (K * p + 1)
  }

  if (spec == "rtrend") {
    LY <- rbind(LY, "trend" = (1 + p):(N - 1))
    LDY <- rbind(LDY, "const" = 1)
    n_par <- K * r + (K + 1) * r + K * (K * p + 1)
  }

  if (spec == "utrend") {
    LDY <- rbind(LDY, "const" = 1, "trend" =  1:n)
    n_par <- K * r + K * r + K * (K * p + 2)
  }
  
  if (!is.null(shift)) {
    n_shift <- length(shift)
    if (n_shift == 1) {
      shift <- c(rep.int(x = 0, times = shift), rep.int(x = 1, times = N - 1 - shift))[-seq_len(p)]
      LY <- rbind(LY, shift)
    } else {
      shift_mat <- matrix(data = NA, nrow = n_shift, ncol = n)
      for (i in seq_len(shift)) {
        shift_mat[i, ] <- c(rep.int(x = 0, times = shift[i]), rep.int(x = 1, times = N - 1 - shift[i]))[-seq_len(p)]
      }
      LY <- rbind(LY, shift_mat)
    }
    n_par <- n_par + r * n_shift
  }

  if (!is.null(season)) {
    if (spec == "none" | spec == "rconst") {
      warning(paste0("'season' is not taken into account because 'spec' = '", spec, "'."), call. = FALSE)
    } else {
      temp <- seasonal_dummies(season = season, season_ref = season_ref)
      season_mat <- temp
      while (nrow(season_mat) < N) {
        season_mat <- rbind(season_mat, temp)
      }
      season_mat <- t(season_mat[(p + 2):N, ])
      LDY <- rbind(LDY, season_mat)
      n_par <- n_par + K * (season - 1)
    }
  }
  
  if (!is.null(exogen_ect)) {
    exogen_ect <- as.matrix(exogen_ect)
    if (NA %in% exogen_ect) {
      stop("'exogen_ect' contains missing values.", call. = FALSE)
    }
    if (nrow(exogen_ect) != N) {
      stop("'exogen_ect' exceeds the number of observations.", call. = FALSE)
    }
    LY <- rbind(LY, t(exogen_ect[-(1:(p + 1)), ]))
    n_par <- n_par + r * ncol(exogen_ect)
  }
  
  if (!is.null(exogen)) {
    exogen <- as.matrix(exogen)
    if (NA %in% exogen) {
      stop("'exogen' contains missing values.", call. = FALSE)
    }
    if (nrow(exogen) != N) {
      stop("'exogen' exceeds the number of observations.", call. = FALSE)
    }
    LDY <- rbind(LDY, t(exogen[-(1:(p + 1)), ]))
    n_par <- n_par + K * ncol(exogen)
  }
  
  for (i in 1:(K + 1)) {
    if (n < n_par[i]) {
      warning("The number of parameters exceeds the number of observations after adjustments.", call. = FALSE)
    }
  }

  n_par <- n_par - r^2

  M <- diag(x = 1, nrow = n, ncol = n) - t(LDY) %*% solve(LDY %*% t(LDY)) %*% LDY

  R_0 <- DY %*% M
  R_1 <- LY %*% M

  S_00 <- R_0 %*% t(R_0) / n
  S_01 <- R_0 %*% t(R_1) / n
  S_10 <- R_1 %*% t(R_0) / n
  S_11 <- R_1 %*% t(R_1) / n

  S_11_sqrt <- sqrtm(S_11)

  S <- solve(S_11_sqrt) %*% S_10 %*% solve(S_00) %*% S_01 %*% solve(S_11_sqrt)

  V <- eigen(x = S, symmetric = TRUE)

  method <- match.arg(method)

  if (method == "johansen") {
    if (is.null(alpha)) {
      alpha <- 0.05
    }
    if (!is.null(shift)) {
      stop("This option is not compatible with shift dummy variables. Please choose 'method' = 'info_criteria'.", call. = FALSE)
    }
    if (!is.null(exogen_ect) || !is.null(exogen)) {
      stop("This option is not compatible with exogenous variables. Please choose 'method' = 'info_criteria'.", call. = FALSE)
    }
    e_val <- V$values[1:K]
    log_e_val <- log(1 - e_val)
    trace_stat <- rep(x = NA, times = K)
    max_stat <- rep(x = NA, times = K)
    for (i in 1:K) {
      trace_stat[i] <- -n * sum(log_e_val[i:K])
      max_stat[i] <- -n * log_e_val[i]
    }
    stat <- cbind("trace_stat" = c(trace_stat, NA), "max_stat" = c(max_stat, NA))
    rownames(stat) <- paste0("r = ", r)
    c_val <- matrix(data = NA, nrow = K + 1, ncol = 2)
    rownames(c_val) <- paste0("r = ", r)
    colnames(c_val) <- c("trace_crit", "max_crit")
    p_val <- matrix(data = NA, nrow = K + 1, ncol = 2)
    rownames(p_val) <- paste0("r = ", r)
    colnames(p_val) <- c("trace_pval", "max_pval")
    for (i in 1:K) {
      output <- johansen_distr(n_sim = 15000, n_step = 1500, rw_dim = K - r[i], spec = spec)
      c_val[i, ] <- apply(X = output, MARGIN = 2, FUN = quantile, probs = 1 - alpha)
      p_val[i, ] <- c(mean(output[, 1] > trace_stat[i]), mean(output[, 2] > max_stat[i]))
    }
    output <- list("n" = n, "p" = p, "spec" = spec, "method" = method, "n_par" = n_par, "e_val" = c(NA, e_val), "stat" = stat, "crit_val" = c_val, "p_val" = p_val, "alpha" = alpha)
    return(structure(.Data = output, class = "cointegration_test"))
  }

  if (method == "info_criteria") {
    det_S_00 <- det(S_00)
    e_val <- V$values[1:K]
    log_e_val <- log(1 - e_val)
    log_sum <- c(0, cumsum(log_e_val))
    log_lik <- rep(x = NA, times = K + 1)
    for (i in 1:(K + 1)) {
      log_lik[i] <- -0.5 * n * (K * (log(2 * pi) + 1) + log(det_S_00) + log_sum[i])
    }
    model_fit <- -2 * log_lik / n
    n_par_res <- n_par / n
    IC <- cbind("AIC" = model_fit + 2 * n_par_res, "BIC" = model_fit + log(n) * n_par_res, "HQC" = model_fit + 2 * log(log(n)) * n_par_res)
    output <- list("n" = n, "p" = p, "spec" = spec, "method" = method, "n_par" = n_par, "e_val" = c(NA, e_val), "log_lik" = log_lik, "IC" = IC)
    return(structure(.Data = output, class = "cointegration_test"))
  }
}
