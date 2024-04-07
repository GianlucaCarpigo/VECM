#' Fit a VEC model
#'
#' This function fits a vector error correction model by using either the maximum likelihood (\code{ML}) or the estimated generalized least squares (\code{EGLS}) method.
#'
#' @param data The matrix containing the multivariate time series to be analyzed.
#' @param p The lag order. See below for more details.
#' @param r The cointegrating rank.
#' @param method The estimation method. It can be selected between \code{ML} (the default) and \code{EGLS}.
#' @param spec The type of deterministic component to be considered. See below for more details.
#' @param shift The break point at which there is a change in the levels of the variables. See below for more details.
#' @param season The number of seasonal dummy variables.
#' @param season_ref The season selected as reference. If \code{NULL}, the last season is automatically used as reference.
#' @param exogen_ect The exogenous variable to be included in the error correction term.
#' @param q_ect The lag order of the exogenous variable in the error correction term.
#' @param exogen The exogenous variable to be included outside the error correction term.
#' @param q The lag order of the exogenous variable outside the error correction term.
#'
#' @details The lag order (\code{p}) to which the user must specify refers to the lag in the VECM representation. It cannot be equal to zero.
#'
#' The seasonal dummies are othogonal and their construction is based on  the description in Lütkepohl(2005).
#'
#' @return An object of class \code{VECM} is a list with the following components: \cr
#'
#' \item{alpha}{The loading matrix.}
#' \item{sigma_alpha}{The covariance matrix of \code{alpha}.}
#' \item{beta}{The cointegrating matrix.}
#' \item{sigma_beta}{The covariance matrix of \code{beta}.}
#' \item{pi}{The long-run impact matrix.}
#' \item{sigma_pi}{The covariance matrix of \code{pi}.}
#' \item{gamma}{The short-run impact matrix.}
#' \item{sigma_gamma}{The covariance matrix of \code{gamma}.}
#' \item{sigma_alphagamma}{The joint covariance matrix of \code{alpha} and \code{gamma}.}
#' \item{sigma_pigamma}{The joint covariance matrix of \code{pi} and \code{gamma}.}
#' \item{u}{The residuals matrix.}
#' \item{sigma_u}{The covariance matrix of \code{u}.}
#' \item{model}{A list containing information about the estimated model.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Johansen, S. (1995). \emph{Likelihood-based inference in cointegrated vector autoregressive models}. Oxford: Oxford University Press.
#'
#' Lütkepohl, H., & Krätzig, M. (2004). \emph{Applied time series econometrics}. Cambridge: Cambridge University Press.
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
#' print(model)
#' summary(model)
#' }
#'
#' @seealso \code{\link{plot.VECM}} \code{\link{print.VECM}} \code{\link{summary.VECM}}
#'
#' @export
#'
#' @importFrom expm sqrtm
#' @importFrom stats embed pt

VECM <- function(data, p, r, method = c("ML", "EGLS"), spec = c("none", "rconst", "uconst", "rtrend", "utrend"), shift = NULL, season = NULL, season_ref = NULL, exogen_ect = NULL, q_ect = NULL, exogen = NULL, q = NULL) {

  data <- as.matrix(data)

  if (NA %in% data) {
    stop("'data' contains missing values.", call. = FALSE)
  }

  K <- ncol(data)
  N <- nrow(data)

  if (is.null(colnames(data))) {
    warning("Variables in 'data' have no name. A predefined labelling is used.", call. = FALSE)
    colnames(data) <- paste0("V", seq_len(K))
  }

  var_names <- colnames(data)

  temp <- seq.int(from = p + 1, to = N - 1)

  LY <- t(data[temp, ])
  rownames(LY) <- paste0(var_names, "_L")

  Z <- t(embed(x = diff(data), dimension = p + 1))

  DY <- Z[seq_len(K), ]
  rownames(DY) <- paste0(var_names, "_D")

  temp <- NULL
  for (i in seq_len(p)) {
    temp <- c(temp, paste0(var_names, "_LD", i))
  }

  LDY <- Z[-seq_len(K), ]
  rownames(LDY) <- temp

  n <- N - (p + 1)

  spec <- match.arg(spec)

  if (spec == "none") {
    D_ect <- NULL
    D <- NULL
    n_par <- K * r + K * r + K^2 * p
  }

  if (spec == "rconst") {
    const <- rep.int(x = 1, times = n)
    LY <- rbind(LY, const)
    D_ect <- rbind(const)
    D <- NULL
    n_par <- K * r + (K + 1) * r + K^2 * p
  }

  if (spec == "uconst") {
    const <- rep.int(x = 1, times = n)
    LDY <- rbind(const, LDY)
    D_ect <- NULL
    D <- rbind(const)
    n_par <- K * r + K * r + K * (K * p + 1)
  }

  if (spec == "rtrend") {
    const <- rep.int(x = 1, times = n)
    trend <- seq.int(from = p + 1, to = N - 1)
    LY <- rbind(LY, trend)
    LDY <- rbind(const, LDY)
    D_ect <- rbind(trend)
    D <- rbind(const)
    n_par <- K * r + (K + 1) * r + K * (K * p + 1)
  }

  if (spec == "utrend") {
    const <- rep.int(x = 1, times = n)
    trend <- seq.int(from = p + 1, to = N - 1)
    LDY <- rbind(const, trend, LDY)
    D_ect <- NULL
    D <- rbind(const, trend)
    n_par <- K * r + K * r + K * (K * p + 2)
  }

  if (!is.null(shift)) {
    n_shift <- length(shift)
    if (n_shift == 1) {
      shift <- c(rep.int(x = 0, times = shift), rep.int(x = 1, times = N - 1 - shift))[-seq_len(p)]
      LY <- rbind(LY, shift)
      D_ect <- rbind(D_ect, shift)
    } else {
      shift_mat <- matrix(data = NA, nrow = n_shift, ncol = n)
      rownames(shift_mat) <- paste0("shift_", seq_len(shift))
      for (i in seq_len(shift)) {
        shift_mat[i, ] <- c(rep.int(x = 0, times = shift[i]), rep.int(x = 1, times = N - 1 - shift[i]))[-seq_len(p)]
      }
      LY <- rbind(LY, shift_mat)
      D_ect <- rbind(D_ect, shift_mat)
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
      rownames(season_mat) <- colnames(temp)
      LDY <- rbind(LDY, season_mat)
      D <- rbind(D, season_mat)
      n_par <- n_par + K * (season - 1)
    }
  }

  if (!is.null(exogen_ect)) {
    exogen_ect <- as.matrix(exogen_ect)
    if (NA %in% exogen_ect) {
      stop("'exogen_ect' contains missing values.", call. = FALSE)
    }
    #exogen_ect_row <- nrow(exogen_ect)
    if (is.null(q_ect)) {
      #if (exogen_ect_row != N) {
      if (nrow(exogen_ect) != N) {
        stop("The number of observations for the exogenous variable must match the number of observations for the endogenous variables.", call. = FALSE)
      }
    }
    if (is.null(colnames(exogen_ect))) {
      warning("Variables in 'exogen_ect' have no name. A predefined labelling is used.", call. = FALSE)
      colnames(exogen_ect) <- paste0("E", 1:ncol(exogen_ect))
    }
    temp <- rownames(LY)
    temp_1 <- colnames(exogen_ect)
    temp_2 <- rownames(D_ect)
    if (is.null(q_ect)) {
      exogen_ect <- t(exogen_ect[-(1:(p + 1)), ])
      rownames(exogen_ect) <- temp_1
      LY <- rbind(LY, exogen_ect)
      rownames(LY) <- c(temp, temp_1)
      D_ect <- rbind(D_ect, exogen_ect)
      rownames(D_ect) <- c(temp_2, temp_1)
    } else {
      if (ncol(exogen_ect) > 1) {
        return()
      }
      temp_1 <- c(temp_1, paste0(temp_1, "_L", 1:q_ect))
      exogen_ect <- embed(x = exogen_ect, dimension = q_ect + 1)
      #w <- exogen_ect_row - n
      w <- nrow(exogen_ect) - n
      if (w > 0) {
        exogen_ect <- t(exogen_ect[-(1:w), ])
        LY <- rbind(LY, exogen_ect)
        rownames(LY) <- c(temp, temp_1)
        D_ect <- rbind(D_ect, exogen_ect)
        rownames(D_ect) <- c(temp_2, temp_1)
      }
      if (w == 0) {
        exogen_ect <- t(exogen_ect)
        LY <- rbind(LY, exogen_ect)
        rownames(LY) <- c(temp, temp_1)
        D_ect <- rbind(D_ect, exogen_ect)
        rownames(D_ect) <- c(temp_2, temp_1)
      }
      if (w < 0) {
        return()
      }
      rownames(exogen_ect) <- temp_1
    }
    n_par <- n_par + r
  }

  if (!is.null(exogen)) {
    exogen <- as.matrix(exogen)
    if (NA %in% exogen) {
      stop("'exogen' contains missing values.", call. = FALSE)
    }
    #exogen_row <- nrow(exogen)
    if (is.null(q)) {
      #if (exogen_row != N) {
      if (nrow(exogen) != N) {
        stop("The number of observations for the exogenous variable must match the number of observations for the endogenous variables.", call. = FALSE)
      }
    }
    if (is.null(colnames(exogen))) {
      warning("Variables in 'exogen' have no name. A predefined labelling is used.", call. = FALSE)
      colnames(exogen) <- paste0("E", 1:ncol(exogen))
    }
    temp <- rownames(LDY)
    temp_1 <- colnames(exogen)
    temp_2 <- rownames(D)
    if (is.null(q)) {
      exogen <- t(exogen[-(1:(p + 1)), ])
      rownames(exogen) <- temp_1
      LDY <- rbind(LDY, exogen)
      rownames(LDY) <- c(temp, temp_1)
      D <- rbind(D, exogen)
      rownames(D) <- c(temp_2, temp_1)
    } else {
      if (ncol(exogen) > 1) {
        return()
      }
      temp_1 <- c(temp_1, paste0(temp_1, "_L", 1:q))
      exogen <- embed(x = exogen, dimension = q + 1)
      #w <- exogen_row - n
      w <- nrow(exogen) - n
      if (w > 0) {
        exogen <- t(exogen[-(1:w), ])
        LDY <- rbind(LDY, exogen)
        rownames(LDY) <- c(temp, temp_1)
        D <- rbind(D, exogen)
        rownames(D) <- c(temp_2, temp_1)
      }
      if (w == 0) {
        exogen <- t(exogen)
        LDY <- rbind(LDY, exogen)
        rownames(LDY) <- c(temp, temp_1)
        D <- rbind(D, exogen)
        rownames(D) <- c(temp_2, temp_1)
      }
      if (w < 0) {
        return()
      }
      rownames(exogen) <- temp_1
    }
    n_par <- n_par + K
  }

  if (n < n_par) {
    stop("The number of parameters exceeds the number of observations after adjustments.", call. = FALSE)
  }

  LY_row <- nrow(LY)
  LY_names <- rownames(LY)

  t_LY <- t(LY)
  t_LDY <- t(LDY)

  M <- diag(x = 1, nrow = n, ncol = n) - t_LDY %*% solve(LDY %*% t_LDY) %*% LDY

  R_0 <- DY %*% M
  R_1 <- LY %*% M
  R_12 <- matrix(data = R_1[-(1:r), ], nrow = LY_row - r, ncol = n, dimnames = list(rownames(R_1)[-(1:r)]))

  S_00 <- tcrossprod(x = R_0, y = R_0) / n
  S_01 <- tcrossprod(x = R_0, y = R_1) / n
  S_10 <- tcrossprod(x = R_1, y = R_0) / n
  S_11 <- tcrossprod(x = R_1, y = R_1) / n

  temp <- solve(sqrtm(S_11))

  S <- temp %*% S_10 %*% solve(S_00) %*% S_01 %*% temp

  method <- match.arg(method)

  if(method == "ML"){

    V <- eigen(x = S, symmetric = TRUE)$vectors[, 1:r]
    beta <- t(t(V) %*% temp)
    beta_last <- matrix(data = (beta %*% solve(beta[1:r, ]))[-(1:r), ], nrow = LY_row - r, ncol = r)
    beta_norm <- rbind(diag(x = 1, nrow = r, ncol = r), beta_last)
    dimnames(beta_norm) <- list(LY_names, LY_names[1:r])
    t_beta_norm <- t(beta_norm)

    alpha <- S_01 %*% beta_norm %*% solve(t_beta_norm %*% S_11 %*% beta_norm)
    colnames(alpha) <- paste0("ECT_", 1:r)

    pi <- alpha %*% t_beta_norm

    gamma <- (DY - pi %*% LY) %*% t_LDY %*% solve(LDY %*% t_LDY)

    u <- DY - pi %*% LY - gamma %*% LDY
    sigma_u <- tcrossprod(x = u, y = u) / n

  }

  if (method == "EGLS"){

    R_11 <- matrix(data = R_1[1:r, ], nrow = r, ncol = n, dimnames = list(rownames(R_1)[1:r]))

    pigamma <- cbind(DY %*% t_LY, DY%*% t_LDY)  %*% solve(cbind(rbind( LY %*% t_LY, LDY %*% t_LY), rbind(LY %*% t_LDY, LDY %*% t_LDY)))
    pi <- pigamma[,1:nrow(LY)]
    alpha <- matrix(data = pi[,1:r], nrow = K, ncol = r)
    rownames(alpha) <- paste0(var_names, "_D")
    colnames(alpha) <- paste0("ECT_", 1:r)
    gamma <- pigamma[,-(1:nrow(LY))]

    u <- DY - pi %*% LY - gamma %*% LDY
    sigma_u <- tcrossprod(x = u, y = u) / n

    beta_last <- solve(t(alpha) %*% solve(sigma_u) %*% alpha) %*% t(alpha) %*% solve(sigma_u) %*% (R_0 - alpha %*% R_11) %*% t(R_12) %*% solve(R_12 %*% t(R_12))
    beta_norm <- rbind(diag(x = 1, nrow = r, ncol = r), t(beta_last))
    dimnames(beta_norm) <- list(LY_names, LY_names[1:r])
    t_beta_norm <- t(beta_norm)
  }

  n_par <- n_par - r^2

  d <- trunc(n_par / K)

  sigma_beta <- n / (n - d) * kronecker(X = solve(tcrossprod(x = R_12, y = R_12)), Y = solve(crossprod(x = alpha, y = solve(sigma_u)) %*% alpha), make.dimnames = TRUE)

  omega_beta <- cbind(rbind(t_beta_norm %*% LY %*% t_LY %*% beta_norm, LDY %*% t_LY %*% beta_norm), rbind(t_beta_norm %*% LY %*% t_LDY, LDY %*% t_LDY))

  sigma_alphagamma <- n / (n - d) * kronecker(X = solve(omega_beta), Y = sigma_u, make.dimnames = TRUE)

  temp <- seq_len(K * r)

  sigma_alpha <- sigma_alphagamma[temp, temp]
  sigma_gamma <- sigma_alphagamma[-temp, -temp]

  omega <- cbind(rbind(LY %*% t_LY, LDY %*% t_LY), rbind(LY %*% t_LDY, LDY %*% t_LDY))

  sigma_pigamma <- n / (n - d) * kronecker(X = solve(omega), Y = sigma_u, make.dimnames = TRUE)

  temp <- seq_len(K^2)

  sigma_pi <- sigma_pigamma[temp, temp]

  model <- list("data" = data, "N" = N, "K" = K, "p" = p, "r" = r, "n" = n, "method" = method, "spec" = spec, "shift" = shift, "season" = season, "season_ref" = season_ref, "exogen_ect" = exogen_ect, "q_ect" = q_ect, "exogen" = exogen, "q" = q, "n_par" = n_par, "d" = d, "DY" = DY, "LY" = LY, "LDY" = LDY, "D_ect" = D_ect, "D" = D, "S" = S, "R_0" = R_0, "R_12" = R_12)

  output <- list("alpha" = alpha, "sigma_alpha" = sigma_alpha, "beta" = beta_norm, "sigma_beta" = sigma_beta, "pi" = pi, "sigma_pi" = sigma_pi, "gamma" = gamma, "sigma_gamma" = sigma_gamma, "sigma_alphagamma" = sigma_alphagamma, "sigma_pigamma" = sigma_pigamma, "u" = u, "sigma_u" = sigma_u, "model" = model)

  return(structure(.Data = output, class = "VECM"))

}
