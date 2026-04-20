#' @title Summary method for 'oscars' objects
#'
#' @description Prints a summary of 'oscars' objects.
#'
#' @inheritParams print.oscars
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
                    , object$evaluations, "evaluations with")
    mess2 <- paste0("functional tolerance stopping rule of ", object$controls$fTol
                  , " and parameter tolerance stopping rule of"
                  , object$controls$xTol)
  }
  mess <- strwrap( c(mess1, mess2) )
  cat(paste(mess, "\n"))
  invisible(NULL)
}
