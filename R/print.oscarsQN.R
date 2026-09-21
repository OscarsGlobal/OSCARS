#' @title Print method for 'oscars' objects
#'
#' @description Prints an 'oscars' object showing minimized (or maximized)
#' parameters and the optimization message.
#'
#' @param x An 'oscars' object returned by \code{oscars}.
#'
#' @param \dots Included for compatibility with other print methods.
#' Ignored here.
#'
#' @return No return value, called for side effects. Technically, \code{NULL}
#' is returned invisibly.
#'
#' @seealso \code{\link{oscars}}
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
#' out <- oscarsQN(hosaki, hosakigrad, 2, 0, upr = c(5,6))
#' out
#'
#' @export
#'
print.oscarsQN <- function(x, ...){

  if( x$controls$DoMax ){
    upDwn <- "Maximum"
  } else {
    upDwn <- "Minimum"
  }
  if( x$convergence == 0 ){
    mess <- paste0(upDwn, " value found at ", paste(x$par, collapse = ", "))
  } else {
    mess <- paste0(upDwn, " not found in ", x$evaluations, " evaluations.")
  }
  mess <- strwrap( mess )
  cat(paste(mess, "\n"))
  invisible(NULL)
}
