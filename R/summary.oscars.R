#' @title Summary method for 'oscars' objects
#'
#' @description Prints a summary of 'oscars' objects.
#'
#' @inheritParams print.oscars
#'
#' @param object An 'oscars' object returned by \code{oscars}.
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
#' summary(out)
#'
#' @export
#'
summary.oscars <- function(object, ...){

  if( x$controls$DoMax ){
    upDwn <- "maximum"
  } else {
    upDwn <- "minimum"
  }
  if( x$convergence == 0 ){
    mess1 <- paste0("Function ", upDwn
                 , " of ", x$value
                 , " achieved in ", x$evaluations, " evaluations"
                 , " at ", paste(x$par, collapse = ", "))
    mess2 <- paste0("Evaluations stopped because ", x$message, ".")
  } else {
    mess1 <- paste0("FAILURE: Function ", upDwn, " not found in "
                    , x$evaluations, "evaluations with")
    mess2 <- paste0("functional tollerance stopping rule of ", x$controls$fTol
                  , " and parameter tollerange stopping rule of"
                  , x$controls$xTol)
  }
  mess <- strwrap( c(mess1, mess2) )
  cat(paste(mess, "\n"))
  invisible(NULL)
}
