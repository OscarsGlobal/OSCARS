#' @title Control parameters for oscars routine
#'
#' @description
#' Provides control over oscar parameters, such as number of iterations,
#' tolerance, etc.
#'
#' @param nfmax The maximum number of function evaluations to perform.
#' Default for \code{nfmax} is 50000.
#'
#' @param infol Verbosity during iterations. If \code{infol} is positive,
#' each new best function value is printed.   Default is 1.
#'
#' @param DoMax logical variable set to TRUE if the objective is
#' to be maximized.   Default is FALSE.
#'
#' @param fTol Stopping tolerance for the objective function f. This tolerance
#' is multiplied by the larger of the absolute value of the current objective
#' function value and 1.  This gives a relative tolerance for large f, and
#' and absolute one otherwise.   If \code{fTol} is negative the algorithm
#' will do the maximum \code{nfmax} of allowed function evaluations and then
#' stop. Default is 1e-6.
#'
#' @param xTol Tolerance in the decision variables which is used to define the
#' minimum sampling box size along each axis.   For each decision variable
#' \code{xTol} is scaled by the larger of 1 and the magnitude of the current
#' value of that variable.   This yields a relative tolerance of \code{xTol}
#' for large magnitude decision variables, and an absolute tolerance for
#' small ones.   Once the sampling box is less than tolerance
#' along all axes, the sequence of nested sample boxes is
#' ended.   \code{xTol} must be positive and the algorithm will impose a
#' minimum value of 1e-12.   Difference between current and previous best
#' known points must be within relative or absolute tolerance of \code{xTol}
#' for oscars to halt before the function budget is exhausted.  Default is 1e-8.
#'
#' @details
#' A subset of parameters can be specified.  All non-specified parameters
#' revert to their defaults.  No parameter abbreviations.
#'
#' @return A named list of control parameters for oscars.
#'
#' @examples
#' oscars.control()  # default values
#' oscars.control(nfmax = 100000) # bump iteration budget
#' oscars.control(xTol = 10*oscars.control()$xTol) # increase xTol
#'
#' @export
oscars.control <- function(
    nfmax = 50000
    , infol = 2
    , DoMax = FALSE
    , fTol  = 1e-6
    , xTol  = 1e-8
){
  contr.list <- list(
    nfmax = nfmax,
    infol = infol,
    DoMax = DoMax,
    fTol = fTol,
    xTol = xTol
  )
  return(contr.list)
}

