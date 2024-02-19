#' Plotting VEC model forecasts
#'
#' This function is a method for the class \code{VECM_forecast}.
#'
#' @param x An object of class \code{VECM_forecast}.
#' @param type The type of plot. It can be chosen between \code{area} (the default), \code{interval}, and \code{line}.
#' @param show The number of observations for which the forecast is computed.
#' @param observed A logical value. If \code{TRUE} (not the default), the observed values are plotted together with the forecasts.
#'
#' @author Gianluca Carpigo \email{gianluca.carpigo@@uniroma1.it}
#'
#' Michele Cianfriglia \email{michele.cianfriglia@@libero.it}
#'
#' @seealso \code{\link{VECM_forecast}}
#'
#' @export
#'
#' @importFrom graphics points polygon title
#' @importFrom utils tail

plot.VECM_forecast <- function(x, type = c("area", "interval", "line"), show = 50, observed = FALSE) {
  if (!inherits(x = x, what = "VECM_forecast")) {
    stop("The object is not of class 'VECM_forecast'.", call. = FALSE)
  }
  data_hat <- rbind(tail(x = x$data - t(x$u), n = show), x$forecast)
  var_names <- colnames(x$data)
  type <- match.arg(type)
  if (type == "area") {
    if (x$type == "static") {
      stop("For a static forecast please select 'type' = 'interval'.", call. = FALSE)
    }
    h <- nrow(x$forecast)
    for (i in 1:length(var_names)) {
      limits <- range(data_hat[, i])
      y_lim <- c(min(limits[1], x$CI_forecast[h, 1, i]) * 0.9, max(limits[2], x$CI_forecast[h, 2, i]) * 1.1)
      x_lim <- tail(x = 1:(nrow(x$data) + h), n = show + h)[c(1, show + h)]
      plot(x = NULL, y = NULL, xlim = x_lim, ylim = y_lim, xlab = "", ylab = "")
      title(xlab = "time", line = 2.3)
      title(ylab = var_names[i], line = 2.3)
      if (observed) {
        points(x = x_lim[1]:nrow(x$data), y = x$data[x_lim[1]:nrow(x$data), i], type = "l", lend = 1, lwd = 2)
      }
      x_1 <- nrow(x$data) + 1
      x_2 <- nrow(x$data) + h
      y_1 <- x$CI_forecast[, 1, i]
      y_2 <- x$CI_forecast[, 2, i]
      polygon(x = c(x_1:x_2, rev(x_1:x_2)), y = c(y_1, rev(y_2)), border = NA, col = "lightgrey")
      points(x = x_lim[1]:x_lim[2], y = data_hat[, i], type = "l", lend = 1, lwd = 2, col = "red")
    }
  }
  if (type == "interval") {
    if (x$type == "static") {
      for (i in 1:length(var_names)) {
        limits <- range(data_hat[, i])
        y_lim <- c(min(limits[1], x$CI_forecast[1, i]) * 0.9, max(limits[2], x$CI_forecast[2, i]) * 1.1)
        x_lim <- tail(x = 1:(nrow(x$data) + 1), n = show + 1)[c(1, show + 1)]
        plot(x = NULL, y = NULL, xlim = x_lim, ylim = y_lim, xlab = "", ylab = "")
        title(xlab = "time", line = 2.3)
        title(ylab = var_names[i], line = 2.3)
        if (observed) {
          points(x = x_lim[1]:nrow(x$data), y = x$data[x_lim[1]:nrow(x$data), i], type = "l", lend = 1, lwd = 2)
        }
        points(x = x_lim[1]:x_lim[2], y = data_hat[, i], type = "l", lend = 1, lwd = 2, col = "red")
        lines(x = c(nrow(x$data) + 1, nrow(x$data) + 1), y = c(x$CI_forecast[1, i], x$CI_forecast[2, i]), type = "l", lend = 1, lwd = 2, col = "orange")
      }
    } else {
      h <- nrow(x$forecast)
      for (i in 1:length(var_names)) {
        limits <- range(data_hat[, i])
        y_lim <- c(min(limits[1], x$CI_forecast[h, 1, i]) * 0.9, max(limits[2], x$CI_forecast[h, 2, i]) * 1.1)
        x_lim <- tail(x = 1:(nrow(x$data) + h), n = show + h)[c(1, show + h)]
        plot(x = NULL, y = NULL, xlim = x_lim, ylim = y_lim, xlab = "", ylab = "")
        title(xlab = "time", line = 2.3)
        title(ylab = var_names[i], line = 2.3)
        if (observed) {
          points(x = x_lim[1]:nrow(x$data), y = x$data[x_lim[1]:nrow(x$data), i], type = "l", lend = 1, lwd = 2)
        }
        points(x = x_lim[1]:x_lim[2], y = data_hat[, i], type = "l", lend = 1, lwd = 2, col = "red")
        for (j in 1:h) {
          lines(x = c(nrow(x$data) + j, nrow(x$data) + j), y = c(x$CI_forecast[j, 1, i], x$CI_forecast[j, 2, i]), type = "l", lend = 1, lwd = 2, col = "orange")
        }
      }
    }
  }
  if (type == "line") {
    if (x$type == "static") {
      stop("For a static forecast please select 'type' = 'interval'.", call. = FALSE)
    }
    h <- nrow(x$forecast)
    for (i in 1:length(var_names)) {
      limits <- range(data_hat[, i])
      y_lim <- c(min(limits[1], x$CI_forecast[h, 1, i]) * 0.9, max(limits[2], x$CI_forecast[h, 2, i]) * 1.1)
      x_lim <- tail(x = 1:(nrow(x$data) + h), n = show + h)[c(1, show + h)]
      plot(x = NULL, y = NULL, xlim = x_lim, ylim = y_lim, xlab = "", ylab = "")
      title(xlab = "time", line = 2.3)
      title(ylab = var_names[i], line = 2.3)
      if (observed) {
        points(x = x_lim[1]:nrow(x$data), y = x$data[x_lim[1]:nrow(x$data), i], type = "l", lend = 1, lwd = 2)
      }
      points(x = x_lim[1]:x_lim[2], y = data_hat[, i], type = "l", lend = 1, lwd = 2, col = "red")
      x_1 <- nrow(x$data) + 1
      x_2 <- nrow(x$data) + h
      points(x = x_1:x_2, y = x$CI_forecast[, 1, i], type = "l", lend = 1, lwd = 2, col = "orange")
      points(x = x_1:x_2, y = x$CI_forecast[, 2, i], type = "l", lend = 1, lwd = 2, col = "orange")
    }
  }
}
