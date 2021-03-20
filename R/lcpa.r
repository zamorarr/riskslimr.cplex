#' Lattice Cutting Plane Algorithm using CPLEX
#'
#' @param x feature matrix
#' @param y response vector
#' @param R_max max number of features (including intercept)
#' @param time_limit max running time in seconds
#' @param logfile file to write log output
#' @export
lcpa_cplex <- function(x, y, logfile, R_max, time_limit) {
  lcpa_cpp(x, y, logfile, R_max, time_limit)
}
