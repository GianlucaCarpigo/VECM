#' Orthogonal complement
#'
#' Internal function used to compute the orthogonal complement of a matrix.
#'
#' @param mat The matrix for which the orthogonal complement is calculated.
#' @param r The rank of the matrix \code{mat}.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @keywords internal

orth_compl <- function(mat, r) {
  mat <- as.matrix(mat)
  n_row <- nrow(mat)
  row_names <- rownames(mat)
  if (r == 0) {
    mat_orth <- diag(x = 1, nrow = n_row, ncol = n_row)
    dimnames(mat_orth) <- list(row_names, row_names)
    return(mat_orth)
  }
  if (r >= 1 & r < n_row) {
    mat_orth <- as.matrix(qr.Q(qr = qr(mat), complete = TRUE)[, -(1:r)])
    dimnames(mat_orth) <- list(row_names, NULL)
    return(mat_orth)
  }
  if (r == n_row) {
    mat_orth <- matrix(data = 0, nrow = n_row, ncol = 1)
    dimnames(mat_orth) <- list(row_names, NULL)
    return(mat_orth)
  }
}
