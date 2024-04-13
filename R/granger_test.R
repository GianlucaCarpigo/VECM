#' Granger-causality test
#'
#' This function implements the Granger-causality test for a VEC model.
#'
#' @param object An object of class \code{VECM}.
#' @param dip A character vector containing the names of the dependent variables.
#' @param indip A character vector containing the names of the independent variables.
#'
#' @return An object of class \code{granger_test} is a list with the following components: \cr
#'
#' \item{R}{The matrix that identifies the parameters to be tested.}
#' \item{stat}{The test statistic.}
#' \item{stat_pval}{The p-value of the test statistic.}
#' \item{df}{The degrees of freedom of the test statistic.}
#' \item{dip}{A character vector containing the names of the dependent variables.}
#' \item{indip}{A character vector containing the names of the independent variables.}
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @references Lütkepohl, H., & Krätzig, M. (2004). \emph{Applied time series econometrics}. Cambridge: Cambridge University Press.
#'
#' Lütkepohl, H. (2005). \emph{New introduction to multiple time series analysis}. New York: Springer.
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
#' test_G <- granger_test(object = model, dip = "R", indip = "Dp")
#' print(test_G)
#'}
#'
#' @export
#'
#' @importFrom stats pchisq

granger_test <- function(object, dip, indip) {

  if (!inherits(object, what = "VECM")) {
    stop("The object is not of class 'VECM'.", call. = FALSE)
  }

  var_names <- colnames(object$model$data)
  
  if (!min(dip %in% var_names) | !min(indip %in% var_names)) {
    stop("The selected variables are not included in the dataset.", call. = FALSE)
  }

  n <- object$model$n
  p <- object$model$p
  pi <- object$pi
  gamma <- object$gamma
  sigma_theta <- object$sigma_pigamma

  temp <- c()
  for (i in 1:ncol(pi)) {
    temp <- c(temp, paste(rownames(pi), "<-", colnames(pi)[i]))
  }
  for (i in 1:ncol(gamma)) {
    temp <- c(temp, paste(rownames(gamma), "<-", colnames(gamma)[i]))
  }

  theta <- matrix(data = cbind(pi, gamma), nrow = nrow(sigma_theta), ncol = 1, dimnames = list(temp))

  link <- NULL
  for (i in 1:length(indip)) {
    link <- c(link, paste0(dip, "_D <- ", indip[i], "_L"))
    for (j in 1:p) {
      link <- c(link, paste0(dip, "_D <- ", indip[i], "_LD", j, sep = ""))
    }
  }

  R <- NULL
  for (i in 1:nrow(theta)) {
    if (rownames(theta)[i] %in% link) {
      R_temp <- rep(x = 0, times = nrow(theta))
      R_temp[i] <- 1
      R <- rbind(R, R_temp)
    }
  }
  rownames(R) <- link

  W_stat <- t(R %*% theta) %*% solve(R %*% sigma_theta %*% t(R)) %*% R %*% theta
  W_df <- length(dip) * length(indip) * (p + 1)
  W_pval <- pchisq(q = W_stat, df = W_df, lower.tail = FALSE)

  output <- list("R" = R, "stat" = W_stat, "stat_pval" = W_pval, "df" = W_df, "dip" = dip, "indip" = indip)

  return(structure(.Data = output, class = "granger_test"))

}
