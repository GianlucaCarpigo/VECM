#' Seasonal dummy variables
#'
#' This function creates seasonal dummy variables.
#'
#' @param season The number of seasons.
#' @param season_ref The season taken as reference.
#'
#' @return A matrix containing the seasonal dummy variables.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @keywords internal

seasonal_dummies <- function(season, season_ref = NULL) {
  mat <- diag(x = 1, nrow = season, ncol = season)
  mat_names <- paste0("season_", seq_len(season))
  if (is.null(season_ref)) {
    mat <- as.matrix((mat - 1 / season)[, -season])
    colnames(mat) <- mat_names[-season]
    return(mat)
  } else {
    mat <- as.matrix((mat - 1 / season)[, -season_ref])
    colnames(mat) <- mat_names[-season_ref]
    return(mat)
  }
}
