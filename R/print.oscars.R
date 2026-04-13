#' @title Print method for 'oscars' objects
#'
#' @description Prints an 'oscars' object.
#'
#' @param x An 'oscars' object returned by \code{oscars}.
#'
#' @param \dots Included for compatibility with other print methods.
#' Ignored here.
#'
#' @return NULL
#'
#' @seealso \code{\link{oscars}}
#'
#' @examples
#' # Branins camel function with global minima of f = -1.0316 at
#' # (0.0898,0.7127) and (0.0898,-0.7127) with four other local minima
#' camel <- function(par) {
#'   x = par[1]
#'   y = par[2]
#'   f = 4*x^2 - 2.1*x^4 + (1/3)*x^6 + x*y + 4*(y^4-y^2)
#'   return(f) }
#' out <- oscars(camel, lwr = c(-5,-5), upr = c(5,5))
#' out
#'
#' @export
#'
print.oscars <- function(x, ...){

  if( x$controls$DoMax ){
    upDwn <- "maximum"
  } else {
    upDwn <- "minimum"
  }
  if( x$convergence == 0 ){
    mess <- paste0("Function ", upDwn, " at ", paste(x$par, collapse = ", "))
  } else {
    mess <- paste0("Function ", upDwn, " not found in ", x$evaluations, "evaluations.")
  }
  mess <- strwrap( mess )
  cat(paste(mess, "\n"))
  invisible(NULL)
}
