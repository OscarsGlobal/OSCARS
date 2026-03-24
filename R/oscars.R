#' @title OSCARS global minimization
#'
#' @description
#' Performs global minimization of a general function subject to box constraints
#' using the OSCARS direct search algorithm (<doi:toChrisPaper>).
#' Oscars performs a fixed number of function evaluations and returns the
#' best known point and the function evaluated at that point.
#'
#' @param fname An R function to be minimized. This function must take a vector
#' of parameter values as its first argument, and return a scalar.  Additional
#' arguments can be supplied via \code{...}
#'
#' @param ... Additional parameters supplied to function \code{fname}.
#'
#' @param lwr A vector of lower bounds for the parameters of \code{fname}.
#' Length is taken to be the number of parameters to minimize over. All bounds
#' must be finite.
#'
#' @param upr A vector of upper bounds for the parameters of \code{fname}. All bounds
#' must be finite.
#'
#' @param nfmax The maximum number of function evaluations to perform.
#'
#' @param infol Verbosity during iterations. If \code{infol} is positive,
#' each new best function value is printed.
#'
#' @return A list containing the best known set of parameters found during
#' \code{nfmax} iterations and the minimized function value at the best known
#' point.
#'
#' @examples
#' # "banana" function with global minimum of zero at (a, a^2)
#' rosenbrock <- function(par, a = 1, b = 100){
#'    (a - par[1])^2 + b*(par[2] - par[1]^2)^2
#' }
#'
#' oscars(rosenbrock, lwr = c(-3,-3), upr = c(3,3)) # min is c(1,1)
#' oscars(rosenbrock, lwr = c(-3,-3), upr = c(3,3), infol = 1)
#'
#' oscars(rosenbrock, a = 2, nfmax = 50000, lwr = c(-10,-10), upr = c(10,10)) # min is c(2,4)
#' oscars(rosenbrock, a = 2, nfmax = 50000, lwr = c(-10,-10), upr = c(10,10), infol = 1)
#'
#'
#' @export
oscars <- function(fname
                  , ...
                  , lwr
                  , upr
                  , nfmax = 1000
                  , infol = -1
                  ){

  # Initialize the method.   n is the number of variables, hmin (> 0) is minimum
  # box size after which the method resets.   Default is 1e-6.
  n = length(lwr)
  if (!exists('hmin')) {
    hmin = 1e-6
  }

  # Define the algorithm parameters.  relcut is a flag.
  # Trent recommends moving these to a 'control' list, see 'nlminb'
  cutposn  = 0.9        # cut position A.  The cut passes through c + (1-A).(x-c)
  cutratio = 1/3        # cutratio. Cuts along all axes > cutratio*longest for x-c
  relcut = TRUE         # relA = 1 reduces A for non-longest axis cuts. (FLAG)
  maxcycle1 = 90        # These two parameters give the max cycle length, which is (90)
  maxcycle2 = 30        # maxcycle1 + cyclenumber * maxcycle2  (30)

  # Initialize flags --- these are logical algorithm variables.
  gogo = TRUE               # Set to zero to halt the program
  NextResetBest = TRUE      # if true next cycle reset is to best point

  # get the function value at the centre of the box.
  edges = upr-lwr
  boxlwr = lwr
  boxupr = upr
  xb = (lwr+upr)/2
  fb = fname(xb, ...)
  xc = xb
  fc = fb

  # Set up the counters
  CycleNr = 1
  nf = 1
  CycleLength = 1
  NewCycle = FALSE
  if (infol > 0){
    cat(sprintf("Iteration %8i.  New best f = %12.6g \n",nf,fb))
  }

  while (gogo){
    oldboxupr = boxupr
    oldboxlwr = boxlwr

    # Find the position of the new test point and its function value
    newx = boxlwr + runif(n)*(boxupr - boxlwr)
    newf = fname(newx, ...)
    nf = nf+1
    CycleLength = CycleLength+1

    # If the new point is better than the control, update the control
    if (newf < fc){  # if 1
      fc = newf
      xc = newx
      boxlwr = lwr
      boxupr = upr
    } else {
      # new point is not better than c, so retain c and shrink the box
      maxstep = max(abs(xc-newx))
      mincutlength = cutratio*maxstep
      for (ii in 1:n) {
        if (abs(xc[ii] - newx[ii]) >= mincutlength) {
          AA = cutposn
          if (relcut & (maxstep > hmin)){
            AA = cutposn*abs(xc[ii] - newx[ii])/maxstep
          }
          # Calculate the cut position and shift the box face inwards
          temp = AA*xc[ii] + (1-AA)*newx[ii]
          if  (xc[ii] > newx[ii]){
            boxlwr[ii] = temp
          } else {
            boxupr[ii] = temp
          }
        } # end of if
      }   # end of for
    }     # end of if 1

    # Check for reset due to max length of cycle reached
    if (CycleLength > maxcycle1 + maxcycle2*CycleNr) {
      NewCycle = TRUE
    }
    if (max(boxupr-boxlwr) < hmin) {
      NewCycle = TRUE
    }

    # Perform the SORC strategy and adjust m counters.
    if (NewCycle) {
        boxlwr = lwr
        boxupr = upr
        CycleNr = CycleNr+1
        CycleLength = 0
        NewCycle = FALSE
        if (NextResetBest)   {
            # Even numbered cycle so initial control point is best point
            xc = xb
            fc = fb
            NextResetBest = FALSE
        } else  {
            # Odd numbered cycle so initial control point is random
            xc = lwr + runif(n)*edges
            fc = fname(xc, ...)
            nf = nf + 1
            NextResetBest = TRUE
        }
    } # end of if


    # Update the best point.
    if (fc < fb) {
      fb = fc
      xb = xc
      if (infol > 0)  cat(sprintf("Iteration %8i.  New best f = %12.6g \n",nf,fb))
    }

    #   Check stopping conditions.
    if (nf >= nfmax){
      gogo = FALSE
    }
  } # end of while

  # Put the final objective and decision variable values into a list and return
  solution <- list(par = xb
                 , value = fb
                 , evaluations = nf
                 , convergence = 0
                 , message = "Maximum iterations reached"
                 )
  return(solution)

}   # end of function.

