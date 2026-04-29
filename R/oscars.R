#' @title OSCARS-II bound constrained global optimization
#' @description
#' Performs black-box optimization of a general function subject to
#' bounds on the unknown parameters using a variant of the OSCARS-II
#' algorithm (Price, Reale and Robertson (2021) <doi.org/10.1007/s10898-020-00928-6>).
#' If all bounds
#' are finite, Oscars acts as a global optimization algorithm.  It has been
#' adapted to handle infinite upper and lower bounds, in which case the method
#' has the characteristics of a local method for nonsmooth problems.  Oscars
#' does not use or assume the existence of derivatives of the objective
#' function.  It is a low overhead method for cheaply evaluated black-box
#' functions.
#' Black-box optimization methods for arbitrary functions do not and cannot
#' provide certificates of optimality if halted after a finite amount of time.
#'
#' Oscars is a stochastic direct search method which uses only function values
#' at selected points.   It generates a finite sequence of nested boxes around
#' a control point, and randomly samples each box once, in turn.   A new set
#' of nested boxes is formed if the current set is exhausted or a point better
#' than the control point is found.   In the latter case the better point
#' replaces the control.   Initially the control point is set to the better
#' of an internal initial point and a user supplied start point (if given).
#'
#' From time to time the control is reset alternately to a random point, or
#' to the best known point.   Each reset marks the end of one cycle and the
#' start of the next.   All even numbered cycles start with control points
#' chosen randomly from the feasible region.   All odd numbered cycles (other
#' than the first) set the control point equal to the best known point.
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
#' @param n The number of parameters with which \code{fname} is minimized.
#'
#' @param lwr A vector of lower bounds for the parameters of \code{fname}.
#' If a single value \code{lwr} is supplied, this value will be used for all
#' lower bounds.  Lower bounds of minus infinity are acceptable.  In order to
#' maximize oscars effectiveness, it is recommended that the gap between
#' upper and lower bounds not be unnecessarily wide.
#'
#'
#' @param upr A vector of upper bounds for the parameters of \code{fname}.
#' If a single value \code{upr} is supplied, this value will be used for all
#' upper bounds.   Upper bounds of infinity are acceptable.   It is suggested
#' that Inf be used for parameters that are unbounded above rather than a
#' very large finite number as this signals the method to operate as a local
#' search rather than attempting to cover all values between the bounds.
#'
#' @param start This is an optional start point for the algorithm.  It allows
#' the user to direct the method to a region the user considers promising.
#' If the start point is infeasible (i.e. violates some bounds) the closest
#' feasible point to it is used.   If a single
#' value is provided, it is used for all dimensions.  Default is null.
#'
#' The algorithm also generates an internal initial point as follows.  When
#' all bounds are finite, this is the centre point of the box.   Otherwise each
#' parameter is started at the average of its bounds when both are finite;
#' if one bound is finite, it uses a feasible value near that bound; otherwise
#' it uses the user supplied start value (if one is given) for that parameter,
#' or zero otherwise.
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
#' # Camel function with global minima of f = -1.0316 at
#' # (0.0898,0.7127) and (0.0898,-0.7127) with four other local minima
#' camel <- function(par) {
#'   x = par[1]
#'   y = par[2]
#'   f = 4*x^2 - 2.1*x^4 + (1/3)*x^6 + x*y + 4*(y^4-y^2)
#'   return(f) }
#' out <- oscars(camel, n = 2, lwr = c(-5,-5), upr = c(5,5))
#'
#' # How to use repeated upper and lower bounds.
#' # Bird function in 2 dimensions.  Global minimum = -106.7645367198
#' bird <- function(par)  {
#'   x1 = par[1];  x2 = par[2]
#'   f = sin(x1)*exp((1-cos(x2))^2) + cos(x2)*exp((1-sin(x1))^2) + (x1-x2)^2
#'   return(f)
#' } # end of bird function
#' out <- oscars(bird, 2, -10, 50)
#'
#' # Hosaki function with global minimum of -2.3458 at (4,2) and one local minimum
#' hosaki <- function(par)  {
#'   x = par[1]
#'   y = par[2]
#'   f = (1 - 8*x + 7*x^2 - (7/3)*x^3 + (1/4)*x^4)*y*y*exp(-y)
#'   return(f) }
#' out <- oscars(hosaki, 2, 0, upr = c(5,6))
#'
#' # The proper way to specify control parameters.
#' out <- oscars(hosaki, 2, lwr = c(0,0), upr = c(5,6),
#'   controls = oscars.control(nfmax = 100000, fTol=10*oscars.control()$fTol))
#'
#' # An example of where the full function evaluation budget is used.
#' out <- oscars(hosaki,2,0,5,controls = oscars.control(nfmax=10000,fTol=-1))
#'
#' # how to pass other values to the objective function
#' # Rosenbrocks "banana" function with global minimum of zero at (a, a^2)
#' rosenbrock <- function(par, a = 1, b = 100){
#'    f = (a - par[1])^2 + b*(par[2] - par[1]^2)^2
#'    return(f)
#' }
#' out <- oscars(rosenbrock, 2, -3, 3, a = 0.5)
#'
#' # Providing a user start point to the algorithm.
#' # Weka_1 function with global minimum of 0 wherever x[1] = -1 and a local
#' # minimum at the origin.  Dimension n is arbitrary with bounds -1 <= x <= 2
#' weka_1 <- function(par)   {
#'   f1 = 1 + sqrt( sum( par^2 ))
#'   f2 = 4*par[1] + 4
#'   f = min(f1,f2)
#'   return(f)
#' }
#' out <- oscars(weka_1, 10, -1, 2, start = 0)
#'
#' # An example of the use of infinite bounds.
#' # Active faces is a nonsmooth function with global minimum = 0 at the
#' # origin.   Standard bounds are 0 <= par <= 5.  Solution is on boundary.
#' activefaces <- function(par)  {
#' f1 = max( log( abs(par) + 1) )
#' f2 = log( abs( sum(par) ) + 1)
#' f = max(f1,f2)
#' return(f)
#' }
#' out <- oscars(activefaces, 10, 0, Inf, start = 4)
#'
#' @export
oscars <- function(fname
                  , n
                  , lwr
                  , upr
                  , ...
                  , start = NULL
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

  # xTol is the tolerance on parameters defining the min sampling box size
  xTol = max(1e-12,xTol)

  # If necessary, fill out the bounds to vectors of length n
  if (length(lwr) < n) {
    temp = lwr[1]
    lwr = rep(temp, times = n)
    if (infol > 0)    print('First lower bound repeated.')
  }
  if (length(upr) < n) {
    temp = upr[1]
    upr = rep(temp, times = n)
    if (infol > 0)    print('First upper bound repeated.')
  }

  # Calculate box edge lengths and check non-negative.
  edges = upr-lwr
  if (min(edges) < 0) {
    if (infol > 0)    print('ALERT: Run aborted due to inconsistent bounds.')
    # Abort run and exit as feasible region is empty.
    solution <- list(par = NULL
                     , value = NULL
                     , evaluations = 0
                     , convergence = 2
                     , message = "ALERT: Run aborted due to inconsistent bounds."
                     , controls = controls
    )
    class(solution) <- "oscars"
    return(solution)
  }

  # Check if the feasible box defined by the bounds is bounded or not.
  BoundedBox = TRUE
  if (max(edges) > 10^300)  BoundedBox = FALSE

  # Define algorithm parameters.  DEFAULT SETTINGS STRONGLY RECOMMENDED.
  # cutposn = cut position A. Cut is through c + (1-A).(x-c) where c is the
  # control point and x is the new point.   0 < A < 1.  Default A = 0.9
  cutposn  = 0.9
  # cutratio. Cuts sampling box along an axis if absolute value of component of
  # x-c along that axis is at least cutratio times largest abs value component.
  # Default = 1/3.
  cutratio = 1/3
  # relcut = 1 reduces A for non-longest axis cuts. Default = TRUE
  relcut = TRUE
  # Each cycle automatically reset after (maxcycle1 + cyclenumber * maxcycle2)
  # function evaluations.   Default is both = 30.
  maxcycle1 = 30
  maxcycle2 = 30
  # The algorithm halts due to negligible progress if the number of cycles
  # exceeds MinNrCycles (default = 8), and the best known function values
  # agree within an absolute (if abs(f) < 1) or relative error of fTol for the
  # last 1-(1/stallratio) of all function evals.   Stallratio default = 3/2.
  stallratio = 3/2
  MinNrCycles = 8
  convergetest = 1
  # When some bounds are infinite the MaxBoxSize defines the largest sampling
  # box in lieu of a bounded feasible region.   (Default = 2)
  MaxBoxSize = 2

  # Initialize logical algorithm variables.
  gogo = TRUE               # Set to zero to halt the program
  NewCycle = FALSE          # Set to true to start a new cycle.
  NextResetBest = FALSE     # if true next cycle reset is to best point
  nanDetected = FALSE       # warns if NaN or NA was found when infol > 0

  # Use all nfmax function evaluations if objective tolerance is negative
  if (fTol < 0) Use_fTol = FALSE  else  Use_fTol = TRUE

  # Check if start point given. If so force feasibility and calculate function.
  if (!is.null(start)) {
    # fill the start point out to a vector if necessary
    if (length(start) < n) {
      temp = start[1]
      start = rep(temp, times = n)
      if (infol > 0)    print('First user start value repeated.')
    }
    # Make sure it satisfies the bounds
    if ((infol > 0) & (min(start-lwr,upr-start) < 0))  {
      print('ALERT: infeasible user supplied start')
    }
    start = pmax(lwr,pmin(upr,start))
    # and then calculate the objective function there
    fstart = fname(start, ...)
    if (DoMax)  fstart = -fstart
    if (is.nan(fstart) | is.na(fstart)) {
      fstart = Inf
      nanDetected = TRUE
    }
  }

  # Form the internal initial point (= centre of box when box is bounded).
  if (BoundedBox) {
    xb = (lwr + upr)/2
  } else {
    # Set up the internal xb for the local search case (not all bounds finite).
    if (!is.null(start)) {
      xb = start
    } else {
      xb = rep(0, times = n)
    }
    # Adjust if one or both bounds are finite.
    for (j in 1:n) {
      if (upr[j] < 10^300)  {
        if (lwr[j] > -10^300)  {
          # Both bounds finite, use average of both bounds.
          xb[j] = (lwr[j]+upr[j])/2
        } else {
          # Unbounded below, step downwards off finite upper bound.
          xb[j] = upr[j] - max(1,abs(upr[j])/2)
        }
      } else {
        if (lwr[j] > -10^300)  {
          # Unbounded above, step upwards off finite lower bound.
          xb[j] = lwr[j] + max(1,abs(lwr[j])/2)
        }
      }
    }
  }
  # Calculate the function at the internal initial point
  fb = fname(xb, ...)
  if (DoMax)  fb = -fb
  if (is.nan(fb) | is.na(fb)) {
    fb = Inf
    nanDetected = TRUE
  }

  # If user start point was given, use if better than internal initial point.
  InternalStart = TRUE
  if (!is.null(start)) {
    if (fstart < fb)  {
      xb = start
      fb = fstart
      if (infol > 0)   cat(sprintf("\n User supplied start:"))
      InternalStart = FALSE
    }
  }
  if (InternalStart & (infol > 0))  cat(sprintf("\n Internal start:"))

  # Now set up the initial control point and second best point
  xc = xb
  fc = fb
  xsecondb = xb
  fsecondb = Inf

  # Set up the mark to monitor progress of f value for fTol stopping rule
  # fmark is restricted to a large finite value to avoid "Inf - Inf" issues.
  fmark = min(10^300,fb)
  nfmark = 1
  fgap = fTol*max(1,abs(fmark))

  # Set up the counters
  CycleNr = 1
  if (is.null(start))  nf = 1   else   nf = 2
  CycleLength = 1
  # Print out function value at better of user and internal initial points
  if (infol > 0){
    if (DoMax) printf = -fb   else   printf = fb
    cat(sprintf("   f = %12.6g \n\n",printf))
  }

  # Set up the bounds on the sampling box.
  if (BoundedBox)  {
    boxlwr = lwr
    boxupr = upr
  } else {
    # box is unbounded so just use large box around current best point xb.
    boxlwr = pmax(lwr,xb - MaxBoxSize*pmax(1,abs(xb)))
    boxupr = pmin(upr,xb + MaxBoxSize*pmax(1,abs(xb)))
  }

  # Begin the main loop.
  while (gogo){
    # Find the position of the new test point and its function value
    newx = boxlwr + stats::runif(n)*(boxupr - boxlwr)
    newf = fname(newx, ...)
    if (DoMax)  newf = -newf
    if (is.nan(newf) | is.na(newf))  {
      newf = Inf
      nanDetected = TRUE
    }
    nf = nf+1
    CycleLength = CycleLength+1

    # Update second best point used for xTol stopping condition.  If newf is
    # new best xb and xsecondb will be corrected when best is updated.
    if (newf < fsecondb) {
      fsecondb = newf
      xsecondb = newx
    }

    if (newf < fc){  # if 1
      # If new point better than control, update control & reset sampling box
      fc = newf
      xc = newx
      if (BoundedBox)  {
        boxlwr = lwr
        boxupr = upr
      } else {
        boxlwr = pmax(lwr,xb - MaxBoxSize*pmax(1,abs(xc)))
        boxupr = pmin(upr,xb + MaxBoxSize*pmax(1,abs(xc)))
      }
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
      # Reset sampling box to feasible box (if bounded), large box otherwise.
      if (BoundedBox)  {
        boxlwr = lwr
        boxupr = upr
      } else {
        boxlwr = pmax(lwr,xb - MaxBoxSize*pmax(1,abs(xb)))
        boxupr = pmin(upr,xb + MaxBoxSize*pmax(1,abs(xb)))
      }
      CycleNr = CycleNr+1
      CycleLength = 0
      NewCycle = FALSE
      if (NextResetBest)   {
        # Odd numbered cycle so initial control point is best point
        xc = xb
        fc = fb
        NextResetBest = FALSE
      } else  {
        # Even numbered cycle.   Initial control point randomly chosen from
        # large sampling box.  (Bounded => lwr = boxlwr and upr = boxupr.)
        xc = boxlwr + stats::runif(n)*(boxupr - boxlwr)
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
      # Print out the new best point information if requested.
      if (infol > 1)  {
        if (DoMax) printf = -fb   else   printf = fb
        cat(sprintf("New best f =   %12.6g   ",printf))
        cat(sprintf("at %7i fevals   in cycle %5i   ",nf,CycleNr))
        cat(sprintf("df = %10.4g   ",abs(fb-fsecondb)))
        dx = sqrt( sum( (xb-xsecondb)^2 ))
        cat(sprintf("dpar = %10.4g   \n",dx))
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
      cat(sprintf("\n Max f = "))
    } else {
      cat(sprintf("\n Min f = "))
    }
    cat(sprintf("%12.6g    Max feval = %7i, used = %7i.   ",fb,nfmax,nf))
    cat(sprintf("Obj Tol = %8.4g   Dec Var Tol = %8.4g \n\n",fTol,xTol))
    if (nanDetected)  {
      cat(sprintf("Note: At least one NaN or NA was returned \n\n"))
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





