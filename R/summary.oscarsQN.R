#' @title Summary method for 'oscarsQN' objects
#'
#' @description Summarizes an 'oscarsQN' object. Shows an 'oscarsQN' object's
#' minimized (or maximized) parameters, optimization message, iterations, etc..
#'
#' @param object An 'oscarsQN' object returned by \code{oscarsQN}.
#'
#' @param ... Ignored here.  Included for use by other methods.
#'
#' @return No return value, called for side effects. Technically, \code{NULL}
#' is returned invisibly.
#'
#' @seealso \code{\link{oscarsQN}}
#'
#' @examples
#' # Hosaki function with global minimum of -2.3458 at (4,2) and one local minimum
#' hosaki <- function(par)  {
#'   x = par[1]
#'   y = par[2]
#'   f = (1 - 8*x + 7*x^2 - (7/3)*x^3 + (1/4)*x^4)*y*y*exp(-y)
#'   return(f) }
#' 
#' hosakigrad <- function(par)  {
#'   x = par[1]
#'   y = par[2]
#'   g = c(0, 0)
#'   g[1] = (-8 + 14*x - 7*x^2 + x^3)*y*y*exp(-y)
#'   g[2] = (1 - 8*x + 7*x^2 - (7/3)*x^3 + (1/4)*x^4)*(2-y)*y*exp(-y)
#'   return(g) }
#' # out <- oscarsQN(hosaki, hosakigrad, 2, 0, upr = c(5,6))


#' summary(out)
#'
#' @export
#'
summary.oscarsQN <- function(object, ...){

  if( object$controls$DoMax ){
    upDwn <- "Maximum"
  } else {
    upDwn <- "Minimum"
  }
  if( object$convergence == 0 ){
    mess1 <- paste0(upDwn
                 , " value of ", object$value
                 , " achieved in ", object$evaluations, " evaluations"
                 , " at ", paste(object$par, collapse = ", "))
    mess2 <- paste0("Evaluations stopped because ", object$message, ".")
  } else {
    mess1 <- paste0("FAILURE: Function ", upDwn, " not found in "
                    , object$evaluations, "evaluations with ")
    mess2 <- paste0(object$controls$kktstopcount, "KKT points with best "
                  , "function value within tolerance of "
                  , object$controls$fTol, "not found and target of "
                  , object$controls$fTarget, " not met.")
  }
  mess <- strwrap( c(mess1, mess2) )
  cat(paste(mess, "\n"))
  invisible(NULL)
}
