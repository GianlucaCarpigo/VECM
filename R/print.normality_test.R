#' Printing multivariate normality tests
#'
#' This function is a method for the class \code{normality_test}.
#'
#' @param x An object of class \code{normality_test}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{normality_test}}
#'
#' @export

print.normality_test <- function(x) {

  if (!inherits(x, what = "normality_test")) {
    stop("The object is not of class VECM_normality.")
  }

  K <- length(x$b1)
  var_names <- names(x$b1)

  skewness_df <- data.frame("skewness" = c(x$b1, NA),
                            "chi-squared" = c(x$b1_chisq, x$b1_joint),
                            "df" = c(rep(x = 1, times = K), K),
                            "p-value" = c(x$b1_pval, x$b1_joint_pval),
                            row.names = c(var_names, "joint"),
                            check.names = FALSE)

  kurtosis_df <- data.frame("kurtosis" = c(x$b2, NA),
                            "chi-squared" = c(x$b2_chisq, x$b2_joint),
                            "df" = c(rep(x = 1, times = K), K),
                            "p-value" = c(x$b2_pval, x$b2_joint_pval),
                            row.names = c(var_names, "joint"),
                            check.names = FALSE)

  JB_df <- data.frame("Jarque-Bera" = c(x$JB, x$JB_joint),
                      "df" = c(rep(x = 2, times = K), 2 * K),
                      "p-value" = c(x$JB_pval, x$JB_joint_pval),
                      row.names = c(var_names, "joint"),
                      check.names = FALSE)

  cat("VEC residual normality tests\n")
  cat(paste("Orthogonalization:", x$type, "\n"))
  cat("Null hypothesis: residuals are multivariate Normal\n\n")
  printCoefmat(skewness_df, P.values = TRUE, has.Pvalue = TRUE, na.print = "", signif.legend = FALSE)
  cat("\n")
  printCoefmat(kurtosis_df, P.values = TRUE, has.Pvalue = TRUE, na.print = "", signif.legend = FALSE)
  cat("\n")
  printCoefmat(JB_df, P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE)

}
