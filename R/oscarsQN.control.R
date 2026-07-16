#' @title Control parameters for oscarsQN routine
#'
#' @description
#' Provides control over oscarsQN parameters, such as number of iterations,
#' tolerance, etc.
#'
#' @param nfmax The maximum number of function evaluations to perform.
#' Default for \code{nfmax} is 50000, but is normally adjusted  by the user.
#'
#' @param infol Verbosity during iterations. If \code{infol} is positive,
#' each new best function value is printed.   Default is 3.
#'
#' @param DoMax logical variable set to TRUE if the objective is
#' to be maximized.   Default is \code{FALSE}.
#'
#' @param fTol Tolerance for comparing f values at KKT points. This tolerance
#' is multiplied by the larger of the absolute value of the current objective
#' function value and 1.  This gives a relative tolerance for large f, and
#' and absolute one otherwise.   
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
#' @param kktTol Tolerance for the Karush-Kuhn-Tucker conditions to determine
#' a potential local minimizer of the problem.   Default = 1e-5.
#' 
#' @param kktstopcount Number of times a KKT point must be identified with
#' objective function value f being within tolerance of the best known value.
#' Maximum accepted variation in f value governed by \code{fTol}.  Default = 3.
#' 
#' @param fTarget Target value for the objective function which, once reached,
#' halts the algorithm immediately.   Default = NULL, which means no target is 
#' set and the method will not halt by this means.
#' 
#' @param CompareGrad If set to TRUE the algorithm with calculate both the
#' analytic gradient and the finite difference gradient estimate.   Both will
#' be printed out for comparison and checking of the analytic gradient code
#' listed under the function \code{gname}.    Execution will proceed with the
#' finite difference gradient.
#'
#' @details
#' A subset of parameters can be specified.  All non-specified parameters
#' revert to their defaults.  No parameter abbreviations.
#'
#' @return A named list of control parameters for oscars.
#'
#' @examples
#' oscarsQN.control()  # default values
#' oscarsQN.control(nfmax = 100000) # bump iteration budget
#' oscarsQN.control(xTol = 10*oscars.control()$xTol) # increase xTol
#'
#' @export
oscarsQN.control <- function(
    nfmax = 50000
    , infol = 3
    , DoMax = FALSE
    , fTol  = 1e-5
    , xTol  = 1e-8
    , kktTol = 1e-5
    , kktstopcount = 3
    , fTarget = NULL
    , CompareGrad = FALSE
){
  contr.list <- list(
    nfmax = nfmax,
    infol = infol,
    DoMax = DoMax,
    fTol = fTol,
    xTol = xTol,
    kktTol = kktTol,
    kktstopcount = kktstopcount,
    fTarget = fTarget,
    CompareGrad = CompareGrad
  )
  return(contr.list)
}

