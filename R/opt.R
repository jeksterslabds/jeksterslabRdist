#' Estimation via Optimization
#'
#' @description Estimates parameters from a given objective function
#'   using optimization.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param FUN Objective function.
#' @param start_values A vector of starting values
#'   corresponding to the number of parameters to be estimated.
#' @param optim Logical.
#'   If `TRUE`,
#'   the `optim` function is used.
#'   If `FALSE`,
#'   the `nlminb` function is used.
#' @param ... Additional arguments to be passed to the optimizer
#'   and arguments to be passed to the objective function.
#' @importFrom stats optim nlminb
#' @return Returns the output of [`optim()`]
#'   or [`nlminb()`].
#' @export
opt <- function(FUN,
                start_values,
                optim = TRUE,
                ...) {
  if (optim) {
    return(
      optim(
        par = start_values,
        fn = FUN,
        ...
      )
    )
  } else {
    return(
      nlminb(
        start = start_values,
        objective = FUN,
        ...
      )
    )
  }
}
