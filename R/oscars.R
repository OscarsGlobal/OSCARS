#' @title OSCARS-II bound constrained global optimization
#' @description
#' Performs black-box global minimization of a general function subject to 
#' bounds on the unknown parameters using a variant of the OSCARS-II 
#' algorithm (https://doi.org/10.1007/s10898-020-00928-6).   Oscars does not
#' use or even assume the existence of derivatives of the objective function.
#' Black-box optimization methods for arbitrary functions do not and cannot
#' provide certificates of optimality if halted after a finite amount of time.
#' Oscars is a low overhead black-box method for cheap to evaluate functions.
#'
#' Oscars is a stochastic direct search method which uses only function values
#' at selected points.   It generates a finite sequence of nested boxes
#' a control point, and randomly samples each box in turn once.   A new set
#' of nested boxes is formed if the current set is exhausted or a point better
#' than the control point is found.   In the latter case the better point
#' replaces the control.   In the first instance, the control point is set at 
#' the centre point of the feasible region given by the lower and upper bounds.  
#' 
#' From time to time the control is reset alternately to a random point, or 
#' to the best known point.   Each reset marks the end of one cycle and the 
#' start of the next.   All even numbered cycles start with control points
#' chosen randomly from the feasible region.   All odd numbered cycles (other
#' than the first) set the control point equal to the best known point.
#' start at the best known point.   Progress in odd numbered cycles suggests
#' the algorithm is refining a minimizer, progress in even cycles suggests
#' it is jumping from one minimizer to another.
#'
#' Oscars either performs a fixed number of function evaluations, or it
#' halts if progress stalls for a significant period of time.  In both cases it
#' returns the best known point and the function evaluated at that point.
#'
#' @param fname An R function to be minimized. This function must take a vector
#' of parameter values as its first argument, and return a scalar.  Additional
#' arguments can be supplied via \code{...}   Missing (NaN and NA) function
#' values are acceptable as they are replaced with Inf when minimizing
#' (or -Inf when maximizing).
#'
#' @param lwr A vector of lower bounds for the parameters of \code{fname}.
#' Length is taken to be the number of parameters to minimize over. All bounds
#' must be finite.
#'
#' @param upr A vector of upper bounds for the parameters of \code{fname}. All
#' bounds must be finite.
#'
#' @param ... Additional parameters supplied to function \code{fname}.
#'
#' @param controls A list of oscar control parameters, such as iteration
#' budget, tolerance, etc. See \code{\link{oscars.control}} for the full list
#' and descriptions.
#'
#' @return A list containing the best known set of parameters found along with
#' function value at the best known parameters.   The number of function 
#' evaluations used, and reason for halting are also given.
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
#'
#' # Bird function in 2 dimensions.  Global minimum = -106.7645367198
#' bird <- function(par)  {
#'   x1 = par[1];  x2 = par[2]
#'   f = sin(x1)*exp((1-cos(x2))^2) + cos(x2)*exp((1-sin(x1))^2) + (x1-x2)^2
#'   return(f)
#' } # end of bird function
#' out <- oscars(bird, c(-10,-10), c(50,50))
#'
#' # Hosaki function with global minimum of -2.3458 at (4,2) and one local minimum
#' hosaki <- function(par)  {
#'   x = par[1]
#'   y = par[2]
#'   f = (1 - 8*x + 7*x^2 - (7/3)*x^3 + (1/4)*x^4)*y*y*exp(-y)
#'   return(f) }
#' out <- oscars(hosaki, lwr = c(0,0), upr = c(5,6))
#'
#' # The proper way to specify control parameters
#' out <- oscars(hosaki, lwr = c(0,0), upr = c(5,6),
#'   controls = oscars.control(nfmax = 100000, fTol=10*oscars.control()$fTol))
#'
#' # Rosenbrocks "banana" function with global minimum of zero at (a, a^2)
#' rosenbrock <- function(par, a = 1, b = 100){
#'    f = (a - par[1])^2 + b*(par[2] - par[1]^2)^2
#'    return(f)
#' }
#'
#' out <- oscars(rosenbrock, lwr = c(-3,-3), upr = c(3,3), a = 0.5)
#'
#' @export
oscars <- function(fname
                  , lwr
                  , upr
                  , ...
                  , controls = oscars.control()
                  ){

  nfmax = controls$nfmax
  infol = controls$infol
  DoMax = controls$DoMax
  fTol = controls$fTol
  xTol = controls$xTol

  if(infol > 0){
    cat("Control parameters:\n")
    print(do.call(cbind, controls))
  }
  
  # Initialize the method.   n is the number of variables.  xTol is the
  # tolerance on decision variables defining the min sampling box size
  n = length(lwr)
  edges = upr-lwr
  xTol = max(1e-12,xTol)

  # Define algorithm parameters.  DEFAULT SETTINGS STRONGLY RECOMMENDED.
  # cutposn = cut position A. Cut is through c + (1-A).(x-c) where c is the
  # control point and x is the new point.   0 < A < 1.  Default 0.9
  cutposn  = 0.9
  # cutratio. Cuts sampling box along an axis if absolute value of component of
  # x-c along that axis is at least cutratio times largest abs value component.
  # Default = 1/3.
  cutratio = 1/3
  # relcut = 1 reduces A for non-longest axis cuts. Default = TRUE
  relcut = TRUE
  # Each cycle automatically reset after (maxcycle1 + cyclenumber * maxcycle2)
  # function evaluations.   Default is 30 + 30*cyclenumber.
  maxcycle1 = 30
  maxcycle2 = 30
  # The algorithm halts due to negligible progress as measured by fTol if the
  # number of cycles exceeds MinNrCycles (default = 8) and the best known
  # function values agree within an absolute (if best known value < 1) or
  # relative error of fTol for the last (stallratio-1) function evaluations.
  stallratio = 3/2
  MinNrCycles = 8
  convergetest = 1

  # Initialize logical algorithm variables.
  gogo = TRUE               # Set to zero to halt the program
  NewCycle = FALSE          # Set to true to start a new cycle.
  NextResetBest = FALSE     # if true next cycle reset is to best point
  nanDetected = FALSE       # warns if NaN or NA found and infol > 0
  # Use all nfmax function evaluations if objective tolerance is negative
  if (fTol < 0) Use_fTol = FALSE  else  Use_fTol = TRUE

  # Set up the bounds on the sampling box.
  boxlwr = lwr
  boxupr = upr
  
  # get initial function value at centre of the box = start point of cycle 1.
  xb = (lwr+upr)/2
  fb = fname(xb, ...)
  if (DoMax)  fb = -fb
  if (is.nan(fb) | is.na(fb)) {
    fb = Inf
    nanDetected = TRUE
  }
  xc = xb
  xsecondb = xb
  fc = fb
  fsecondb = Inf
  
  # Set up the mark to monitor progress of f value for fTol stopping rule
  # fmark is restricted to a large finite value to avoid "Inf - Inf" issues.
  fmark = min(10^300,fb)
  nfmark = 1
  fgap = fTol*max(1,abs(fmark))

  # Set up the counters
  CycleNr = 1
  nf = 1
  CycleLength = 1
  # Print out the function value at the initial point = centre of box
  if (infol > 0){
    if (DoMax) printf = -fb   else   printf = fb
    cat(sprintf("Iteration %8i.  New best f = %12.6g \n",nf,printf))
  }

  while (gogo){
    # Find the position of the new test point and its function value
    newx = boxlwr + runif(n)*(boxupr - boxlwr)
    newf = fname(newx, ...)
    if (DoMax)  newf = -newf
    if (is.nan(fb) | is.na(fb))  {
      newf = Inf
      nanDetected = TRUE
    }
    nf = nf+1
    CycleLength = CycleLength+1
    
    # Update second best used for xTol stopping condition.  If newf is new
    # best xb and xsecondb will be corrected when best is updated.
    if (newf < fsecondb) {
      fsecondb = newf
      xsecondb = newx
    }

    if (newf < fc){  # if 1
      # If new point better than control, update control & reset sampling box
      fc = newf
      xc = newx
      boxlwr = lwr
      boxupr = upr
    } else {
      # new point not better than control so keep control & shrink sampling box
      maxstep = max(abs(xc-newx))
      mincutlength = cutratio*maxstep
      # shorten the sampling box along each axis for which the control and
      # new point are at least mincutlength apart.
      for (ii in 1:n) {
        if (abs(xc[ii] - newx[ii]) >= mincutlength) {
          AA = cutposn
          if (relcut & (maxstep > xTol)){
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

    # Check for reset due to maximum length of cycle reached
    if (CycleLength > maxcycle1 + maxcycle2*CycleNr) {
      NewCycle = TRUE
    }
    
    # Check for reset if sampling box too small
    abs_rel_xTol = pmax(1,abs(xc))
    if (max(boxupr-boxlwr - xTol*abs_rel_xTol) <= 0) {
      NewCycle = TRUE
    }

    # Start a new cycle.
    if (NewCycle) {
        # Reset sampling box to full feasible region.
        boxlwr = lwr
        boxupr = upr
        CycleNr = CycleNr+1
        CycleLength = 0
        NewCycle = FALSE
        if (NextResetBest)   {
            # Odd numbered cycle so initial control point is best point
            xc = xb
            fc = fb
            NextResetBest = FALSE
        } else  {
            # Even numbered cycle so initial control point is random
            xc = lwr + runif(n)*edges
            fc = fname(xc, ...)
            if (DoMax)  fc = -fc
            if (is.nan(fc) | is.na(fc))  {
              fc = Inf
              nanDetected = TRUE
            }
            nf = nf + 1
            NextResetBest = TRUE
        }
    } # end of if

    # Update the best point.
    if (fc < fb) {
      fsecondb = fb
      fb = fc
      xsecondb = xb
      xb = xc
      if (infol > 0)  {
        if (DoMax) printf = -fb   else   printf = fb
        cat(sprintf("New best f =   %12.6g   ",printf))
        cat(sprintf("at %7i fevals   in cycle %5i   \n",nf,CycleNr))
        }
      # Check if significant progress has been made.  If so reset fmark
      if (Use_fTol & (fb <= fmark - fgap))   {
        fmark = fb
        nfmark = nf
        # update f value gap within absolute or relative fTol of fmark
        fgap = fTol*max(1,abs(fmark))
      }
    }

    #   Check stopping conditions.
    if (nf >= nfmax){
      gogo = FALSE
      message = "Maximum iterations reached"
      convergetest = 0
    }
    if (Use_fTol) {
      if ((fb > fmark - fgap) & (nf > stallratio*nfmark) & (CycleNr > MinNrCycles))  {
        xdiff = abs(xsecondb - xb)
        if (max(xdiff - xTol*pmax(1,abs(xb))) <= 0) {
          gogo = FALSE
          message = "Optimum function value tolerance reached"
          convergetest = 0
        }
      }
    }
  } # end of while

  if (DoMax)  fb = -fb
  if (infol > 0) {
    if (DoMax) { 
      cat(sprintf("\n Max ")) 
    } else { 
      cat(sprintf("\n Min ")) 
    }
    cat(sprintf("problem.  Max feval = %7i, used = %7i.   ",nfmax,nf))
    cat(sprintf("Obj Tol = %8.4g   Dec Var Tol = %8.4g \n\n",fTol,xTol))
    if (nanDetected)  {
      cat(sprintf("Warning: at least one NaN or NA was returned \n\n"))
    }
  }

  # Put the final objective and decision variable values into a list and return
  solution <- list(par = xb
                 , value = fb
                 , evaluations = nf
                 , convergence = convergetest
                 , message = message
                 , controls = controls
  )
  class(solution) <- "oscars"
  return(solution)

}   # end of function.



