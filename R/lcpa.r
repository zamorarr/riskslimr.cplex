#' Lattice Cutting Plane Algorithm using CPLEX
#'
#' @param x feature matrix
#' @param y response vector
#' @param weights vector of observation weights
#' @param R_max max number of features (including intercept)
#' @param time_limit max running time in seconds
#' @param logfile file to write log output
#' @export
lcpa_cplex <- function(x, y, weights, logfile, R_max, time_limit) {
  n <- nrow(x)
  d <- ncol(x)
  stopifnot(identical(length(y), n))
  stopifnot(identical(length(weights), n))

  lcpa_cpp(x, y, weights, logfile, R_max, time_limit)
}
