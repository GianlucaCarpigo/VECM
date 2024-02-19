#' ECT_plot
#'
#' Auxiliary function used in \code{\link{plot.VECM}}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @keywords internal

ECT_plot <- function() {
  cat("Which representation of the error correction terms do you want to plot?\n\n")
  cat("1: representation without deterministic terms\n")
  cat("2: representation with deterministic terms (if any)\n")
  cat("3: representation without deterministic terms and without short-run dynamics\n")
  cat("4: representation with deterministic terms (if any) but without short-run dynamics\n\n")
  cat("Enter your selection:")
  readline()
}

#' VECM_plot
#'
#' Auxiliary function used in \code{\link{plot.VECM}}.
#'
#' @keywords internal

VECM_plot <- function() {
  cat("Make a plot selection or press ESC to exit. Options available:\n\n")
  cat(" 1: error correction terms\n")
  cat(" 2: raw residuals\n")
  cat(" 3: raw residuals standardized\n")
  cat(" 4: squared residuals\n")
  cat(" 5: squared residuals standardized\n")
  cat(" 6: autocorrelation functions of raw residuals\n")
  cat(" 7: autocorrelation functions of squared residuals\n")
  cat(" 8: partial autocorrelation functions of raw residuals\n")
  cat(" 9: partial autocorrelation functions of squared residuals\n")
  cat("10: kernel density estimations of raw residuals\n\n")
  cat("Enter your selection:")
  readline()
}
