#' @title OSCARS-II-quasi-Newton for bound constrained global optimization
#' @description
#' Performs global optimization of a general function subject to
#' bounds on the unknown parameters using a variant of the algorithm in
#' (Price, Reale and Robertson (2025) <doi.org/10.1007/s43069-024-00403-y>).
#' The method is designed for functions which are continuously differentiable,
#' but will run on black-box functions where only function values are available.
#' It incorporates aspects of the direct search method OSCARS-II, guaranteeing 
#' eventual convergence even on black-box continuous functions.
#' From time to time quasi-Newton steps are performed to accelerate convergence
#' and improve the accuracy of the estimated optimizer.   The algorithm will
#' use analytic gradients if provided, otherwise it will estimate them via
#' finite differences.
#' 
#' If all bounds are finite, oscarsQN acts as a global optimization algorithm.  
#' It has been adapted to handle infinite upper and lower bounds, in which case 
#' the method has the characteristics of a local method.  
#' Black-box methods for global optimization of arbitrary functions do not and 
#' cannot provide certificates of optimality if halted after a finite amount 
#' of time, even if the gradient is available at sample points.
#'
#' OscarsQN is a stochastic direct search method which uses function values
#' and gradients at selected points.   It generates a finite sequence of nested 
#' boxes around a control point, and randomly samples each box once, in turn. 
#' Once the current set is exhausted or a point better than the control point 
#' is found the algorithm performs one quasi-Newton step and constructs a new 
#' set of nested boxes.   If a better point than the control is found, it
#' replaces the control.   Initially the control point is set to the better
#' of an internal initial point and a user supplied start point (if given).
#'
#' From time to time the control is reset alternately to a random point, or
#' to the best known point.   Each reset marks the end of one cycle and the
#' start of the next.   
#' 
#' OscarsQN either performs a fixed number of function evaluations, or it
#' halts if a user specified target value has been reached, or the same best
#' locally optimal point has been seen a prescribed number of times.   It
#' returns the best known point and the function value at that point.
#'
#' 
#' @param fname An R function to be minimized. This function must take a vector
#' of parameter values as its first argument, and return a scalar.  Additional
#' arguments can be supplied via \code{...}   Missing (NaN and NA) function
#' values are acceptable as they are replaced with Inf when minimizing
#' (or -Inf when maximizing).
#' 
#' @param gname An R function which returns the gradient vector of the function
#' \code{fname}.  Default is NULL which results in finite differences being
#' used to estimate gradients.   A setting is available in the controls which
#' calculates both the finite difference and analytic gradients for comparison.
#' When both are calculated, the finite difference gradients are used.
#' Additional arguments can be supplied via \code{...}
#'
#' @param n The number of parameters with which \code{fname} is minimized.
#'
#' @param lwr A vector of lower bounds for the parameters of \code{fname}.
#' If a single value \code{lwr} is supplied, this value will be used for all
#' lower bounds.  Lower bounds of minus infinity are acceptable.  In order to
#' maximize oscarsQN's effectiveness, it is recommended that the gap between
#' upper and lower bounds not be unnecessarily wide.  It is suggested
#' that -Inf be used for parameters that are unbounded below rather than a
#' very large (negative) number as this signals the method to act as a local
#' search rather than attempting to cover all values between the bounds.
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
#' @param ... Additional parameters supplied to the functions \code{fname}
#' and \code{gname}.
#'
#' @param controls A list of oscar control parameters, such as iteration
#' budget, tolerance, etc. See \code{\link{oscarsQN.control}} for the full list
#' and descriptions.
#'
#' @return A list containing results of the optimization.  This list consists 
#' of the following components: 
#'   \itemize{
#'      \item \code{par}: vector containing the best known parameters
#'      \item \code{value}: The minimized (or maximized) function value
#'      \item \code{evaluations}: The number of function evaluations used.  
#'            Function evaluations used in calculating finite difference 
#'            estimates of the gradient are included in this total, but
#'            any analytic gradient calculations are not.
#'      \item \code{cycles}: The number of cycles used.
#'      \item \code{convergence}: 0 if the function value target has been
#'         reached, or the required number of best KKT points have been
#'         found; 1 if neither of these cases have not been achieved but 
#'         function evaluation budget has been exhausted; 2 if bounds are 
#'         inconsistent.
#'      \item \code{message}: A text string explaining the value 
#'                  in \code{convergence}.
#'      \item \code{numberKKTpoints}: Number of Karush-Kuhn-Tucker (KKT) 
#'                  points found.   KKT points are potential minimizers.
#'      \item \code{numberbestKKTpoints}: Number of KKT points found which 
#'         take the best known function value (within tolerance).
#'      \item \code{controls}: The values of the controls provided to oscarsQN.
#'      
#' @examples
#' # Camel function with global minima of f = -1.0316 at
#' # (0.0898,0.7127) and (0.0898,-0.7127) with four other local minima
#' camel <- function(par) {
#'   x = par[1]
#'   y = par[2]
#'   f = 4*x^2 - 2.1*x^4 + (1/3)*x^6 + x*y + 4*(y^4-y^2)
#'   return(f) }
#' 
#' camelgrad <- function(par) {
#'   x = par[1]
#'   y = par[2]
#'   g = c(0, 0)
#'   g[1] = 8*x - 8.4*x^3 + 2*x^5 + y
#'   g[2] = x + 16*y^3 - 8*y
#'   return(g) }
#' # out <- oscarsQN(camel, camelgrad, n = 2, lwr = c(-5,-5), upr = c(5,5))
#' 
#' 
#' # Bird function in 2 dimensions.  Global minimum = -106.7645367198
#' bird <- function(par)  {
#'   x1 = par[1];  x2 = par[2]
#'   f = sin(x1)*exp((1-cos(x2))^2) + cos(x2)*exp((1-sin(x1))^2) + (x1-x2)^2
#'   return(f) } 
#' 
#' birdgrad <- function(par)  {
#'   x1 = par[1];  x2 = par[2]
#'   g = c(0, 0)
#'   g[1] = cos(x1)*exp((1-cos(x2))^2) - 2*cos(x2)*exp((1-sin(x1))^2)*(1-sin(x1))*cos(x1) + 2*(x1-x2)
#'   g[2] = 2*sin(x1)*exp((1-cos(x2))^2)*(1-cos(x2))*sin(x2) - sin(x2)*exp((1-sin(x1))^2) + 2*(x2-x1)
#'   return(g) } 
#' # out <- oscarsQN(bird, 2, -10, 50)
#' 
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
#' 
#' # Rosenbrocks "banana" function with global minimum of zero at (a, a^2)
#' rosenbrock <- function(par, a = 1, b = 100) {
#'   f = (a - par[1])^2 + b*(par[2] - par[1]^2)^2
#'   return(f)  }
#' 
#' rosenbrockgrad <- function(par, a = 1, b = 100) {
#'   g = c(0, 0)
#'   g[1] = -2*(a - par[1]) + 2*b*(par[2] - par[1]^2)*(-2*par[1])
#'   g[2] = 2*b*(par[2] - par[1]^2)
#'   return(g)  }
#' # out <- oscarsQN(rosenbrock, rosenbrockgrad, 2, -3, 3, a = 0.5)
#' 
#' # Schwefel function with global min of -418.9829n at x_i = 420.97...
#' # in n dimensions, where n is arbitrary.
#' schwefel <- function(par) {
#'   f = - sum(par*sin(sqrt(abs(par))))
#'   return(f)  }
#' 
#' schwefelgrad <- function(par) {
#'   rootpar = sqrt(abs(par))
#'   g = -sin(rootpar) - 0.5*rootpar*cos(rootpar)
#'   return(g)  }
#' # out <- oscarsQN(schwefel, schwefelgrad, dimension_of_problem, -500, 500)
#' 
#' # vardim function with global min of 0 at par[i] = 1 in n dimensions.
#' vardim <- function(par) {
#'   n = length(par)
#'   temp = c(1:n)
#'   fn1 = sum(temp*(par-1))
#'   f = sum((par-1)^2) + fn1^2 + fn1^4
#'   return(f)   }
#' 
#' vardimgrad <- function(par) {
#'   n = length(par)
#'   temp = c(1:n)
#'   fn1 = sum(temp*(par-1))
#'   g = 2*(par-1) + (2*fn1 + 4*fn1^3)*temp
#'   return(g)  }
#' # out <- oscarsQN(vardim, vardimgrad, dimension_of_problem, 0, 2.7182818)

#' # dixon function with global min of 0 in n dimensions at par[i] = 1.
#' dixon <- function(par) {
#'   n = length(par)
#'   xlo = par[1:n-1]
#'   xhi = par[2:n]
#'   f = (1-par[1])^2 + (1-par[n])^2 + sum( (xlo^2 - xhi)^2 )
#'   return(f)   }
#' 
#' dixongrad <- function(par) {
#'   n = length(par)
#'   x = par
#'   g = rep(0, times = n)
#'   g[1] = -2*(1-x[1]) + 2*(x[1]^2 - x[2])*2*x[1]
#'   jk = 2
#'   while (jk < n)  {
#'     #cat(sprintf("j = %i \n",jk))
#'     g[jk] = 2*(x[jk]^2 - x[jk+1])*2*x[jk] - 2*(x[jk-1]^2 - x[jk])
#'     jk = jk+1
#'   }
#'   g[n] = -2*(1-x[n]) + 2*(x[n-1]^2 - x[n])*(-1)
#'   return(g)  }
#' # out <- oscarsQN(dixon, dixongrad, dimension_of_problem, -2, 2)
#' 
#' @export

oscarsQN <- function(fname
                     , gname = NULL
                     , n
                     , lwr
                     , upr
                     , ...
                     , start = NULL
                     , controls = oscarsQN.control()
){
  
  nfmax = controls$nfmax
  infol = controls$infol
  DoMax = controls$DoMax
  fTol = controls$fTol
  xTol = controls$xTol
  kktTol = controls$kktTol
  KKTstopcount = controls$kktstopcount
  fTarget = controls$fTarget
  CompareGrad = controls$CompareGrad
  
  # Set the tolerance for the Quadratic Programming subproblem stopping rule.
  QPGradTol = 1e-12
  
  # fix the random seed for testing purposes
  set.seed(123)
  
  if(infol > 0){
    cat("Control parameters:\n")
    print(do.call(cbind, controls))
  }

  # xTol is the tolerance on parameters defining the min sampling box size
  xTol = max(1e-12,xTol)
  fTol = max(0,fTol)
  
  # If a target value for f has been given, adjust it for max / min option
  if (is.null(fTarget)) {
    fTarget = -Inf
  } else {
    if (DoMax)   fTarget = -fTarget
  }
  
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
                     , cycles = 0
                     , convergence = 2
                     , message = "ALERT: Run aborted due to inconsistent bounds."
                     , numberKKTpoints = 0
                     , numberbestKKTpoints = 0
                     , controls = controls
    )
    class(solution) <- "oscarsgrad"
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
  
  # Set up the counters
  CycleNr = 1
  if (is.null(start))  nf = 1   else   nf = 2
  CycleLength = 1
  KKTPointCount = 0
  KKTbestcount = 0

  # Set the finite difference increment value
  h = 1e-6
  
  # Now set up the initial control point, and calculate the gradient there.
  xc = xb
  fc = fb
  gradlist <- GetGradient(fname,gname,fc,xc,h,lwr,upr,n,CompareGrad, ...)
  gc = gradlist$gc
  fcount = gradlist$fcount
  nf = nf+fcount
  if (DoMax)  gc = -gc
  gb = gc
  oldxc = xc
  oldgc = gc
  
  # Set up the final point to be returned and initial best kkt point value
  ffinal = fc
  xfinal = xc
  fbestkkt = Inf
  
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
  
  # Set up the initial and maximum Trust region radii
  maxTRR = min( sqrt(sum(edges^2)) , 1e300)
  TRR = maxTRR/2
  TRRbest = TRR
  minTRR = 1e-12
  
  # Set up the initial Hessian estimate as the identity
  B = diag(n)
  Bbest = B
  
  # Beginning of while.   ( BEGINNING OF MAIN LOOP )  ######################
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

    KKTPointFound = FALSE
    NewPass = FALSE      
    QNStepTaken = FALSE
    if (newf < fc){  # if 1
      # New point better than control, so update control and flag end of pass
      oldfc = fc
      oldgc = gc
      oldxc = xc
      fc = newf
      xc = newx
      # calculate the gradient at the new xc
      gradlist <- GetGradient(fname,gname,fc,xc,h,lwr,upr,n,CompareGrad, ...)
      gc = gradlist$gc
      fcount = gradlist$fcount
      nf = nf+fcount
      if (DoMax)  gc = -gc
      NewPass = TRUE
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
    if (max(boxupr-boxlwr - xTol*abs_rel_xTol) <= 0)   NewCycle = TRUE
    
    # Check the current gradient is finite
    gc.ok = CheckGradient(gc,n)
    if (!gc.ok)   nanDetected = TRUE
    
    # Do 1 iteration of Quasi Newton at end of pass or cycle if gc is finite. 
    # if ( (NewCycle) & gc.ok )  {
    if ( (NewPass | NewCycle) & gc.ok )  {
        # First do the BFGS update on B if oldgc is finite.
      oldgc.ok = CheckGradient(oldgc,n)
      if (!oldgc.ok)   nanDetected = TRUE
      if (oldgc.ok)   B <- UpdateB(B,xc,gc,oldxc,oldgc,y,n)

      # Then solve the local box constrained QP for proposed step
      QPsol.out <- CLboxQPsolver(xc,gc,B,TRR,upr,lwr,n,QPGradTol)
      xnew = QPsol.out$xnew
      qnew = QPsol.out$qnew

      # Project new point into feasible region and get its function value 
      xnew = pmin(upr,pmax(lwr,xnew))
      fnew = fname(xnew, ...)
      if (DoMax)  fnew = -fnew
      if (is.nan(fnew) | is.na(fnew))  {
        fnew = Inf
        nanDetected = TRUE
      }
      nf = nf+1

      # Next update control if descent obtained
      if (fnew < fc)    {
        QNStepTaken = TRUE
        # New point better than control, so update control
        oldfc = fc
        oldgc = gc
        oldxc = xc
        fc = fnew
        xc = xnew
        # calculate the gradient at the new xc
        gradlist <- GetGradient(fname,gname,fc,xc,h,lwr,upr,n,CompareGrad, ...)
        gc = gradlist$gc
        fcount = gradlist$fcount
        nf = nf+fcount
        if (DoMax)  gc = -gc
      }

      # Calculate the descent ratio
      DescentRatio = 0
      if (QNStepTaken & qnew < 0)  DescentRatio = (oldfc-fc)/(-qnew)
      
      # Update trust region radius (typical is 0.1 and 0.5, divide by 4 times 2)
      if (DescentRatio < 0.1) {
        TRR = max(minTRR,TRR/4)
      } else if (DescentRatio > 0.5)   {
        TRR = min(maxTRR,2*TRR);
      }
    }
    
    # Update the best point.
    if (fc < fb) {
      fb = fc
      xb = xc
      gb = gc
      Bbest = B
      TRRbest = TRR
      # Update the final point to be reported
      if (fb < ffinal)   {
        ffinal = fb
        xfinal = xb
      }
      # If new control point is a KKT point reset working best to fb = Inf,
      # and start a new cycle from a random control point
      if (CheckGradient(gb,n))   {
        if ( KKTtest(xb,gb,lwr,upr,kktTol) )  {
          KKTPointFound = TRUE
          KKTPointCount = KKTPointCount + 1
          # Update the count of number of times this kkt fbest has been found
          if (fb < fbestkkt - fTol*max(1,abs(fb)))   {
            fbestkkt = fb
            KKTbestcount = 1
          } else if (fb < fbestkkt + fTol*max(1,abs(fb))) {
            fbestkkt = min(fbestkkt,fb)
            KKTbestcount = KKTbestcount + 1
          } 
          # Now reset to the start of a new cycle and reset the working best.
          NewCycle = TRUE
          NextResetBest = FALSE
          fb = Inf
        }
      }
      # Print out the new best point information if requested.
      if (infol > 1)  {
        if (DoMax) printff = -ffinal   else   printff = ffinal
        if (DoMax) printfb = -fb       else   printfb = fb
        cat(sprintf("Best f = %12.6g  ",printff))
        cat(sprintf("Working f = %10.4g  ",printfb))
        cat(sprintf("at %7i fevals, cycle %4i  ",nf,CycleNr))
        cat(sprintf("TRR = %10.3g  ",TRR))
        if (QNStepTaken)     cat(sprintf("q-N step taken."))
        if (KKTPointFound)   cat(sprintf("\n\n KKT point found. \n"))
        cat(sprintf("\n"))
      }
    }
    
    # If new cycle or new pass, reset the sampling box
    if (NewPass | NewCycle)  {
      # Reset sampling box to feasible box (if bounded), large box otherwise.
      if (BoundedBox)  {
        boxlwr = lwr
        boxupr = upr
      } else {
        boxlwr = pmax(lwr,xb - MaxBoxSize*pmax(1,abs(xb)))
        boxupr = pmin(upr,xb + MaxBoxSize*pmax(1,abs(xb)))
      }
    }
    
    # Start a new cycle.
    if (NewCycle) {
      CycleNr = CycleNr+1
      CycleLength = 0
      NewCycle = FALSE
      if (NextResetBest)   {
        # Initial control point of next cycle is working is best point
        xc = xb
        fc = fb
        gc = gb
        B = Bbest
        TRR = TRRbest
        NextResetBest = FALSE
      } else  {
        # Initial control point randomly chosen from
        # large sampling box.  (Bounded => lwr = boxlwr and upr = boxupr.)
        xc = boxlwr + stats::runif(n)*(boxupr - boxlwr)
        fc = fname(xc, ...)
        if (DoMax)  fc = -fc
        if (is.nan(fc) | is.na(fc))  {
          fc = Inf
          nanDetected = TRUE
        }
        nf = nf + 1

        gradlist <- GetGradient(fname,gname,fc,xc,h,lwr,upr,n,CompareGrad, ...)
        gc = gradlist$gc
        fcount = gradlist$fcount
        nf = nf+fcount
        if (DoMax)  gc = -gc
        B = diag(n)
        NextResetBest = TRUE
      }
    } # end of if
    
    #   Check stopping conditions.
    if (ffinal < fTarget)  {
      gogo = FALSE
      message = "Target function value reached."
      convergetest = 0
    }
      if (KKTbestcount >= KKTstopcount) {
      gogo = FALSE
      message = "Best value optimal within tolerance required number of times"
      convergetest = 0
    }
    if (nf >= nfmax){
      gogo = FALSE
      message = "Maximum iterations reached"
      convergetest = 1
    }
  } # end of while  (END OF MAIN LOOP)   #####################################
  
  if (DoMax)  ffinal = -ffinal
  if (infol > 0) {
    if (DoMax) {
      cat(sprintf("\n Max f = "))
    } else {
      cat(sprintf("\n Min f = "))
    }
    cat(sprintf("%12.6g    Max feval = %7i, used = %7i.   ",ffinal,nfmax,nf))
    cat(sprintf("KKT Tol on f = %8.4g   Dec Var Tol = %8.4g \n\n",fTol,xTol))
    cat(sprintf(" Number of KKT points found = %2i   Number ",KKTPointCount))
    cat(sprintf("within tolerance of best known = %2i \n\n",KKTbestcount))
    if (nanDetected)  {
      cat(sprintf("Note: At least one NaN or NA was returned \n\n"))
    }
  }
  
  # Put the final objective and decision variable values into a list and return
  solution <- list(par = xfinal
                   , value = ffinal
                   , evaluations = nf
                   , cycles = CycleNr
                   , convergence = convergetest
                   , message = message
                   , numberKKTpoints = KKTPointCount
                   , numberbestKKTpoints = KKTbestcount
                   , controls = controls
  )
  class(solution) <- "oscars"
  return(solution)
  
}   # end of function.

#------------------------------------------------------------------------------
KKTtest <- function(xx,gg,lwr,upr,kktTol)    {
# calculates the length of the projected reduced gradient.   If below 
# tolerance the KKT conditions hold within tolerance

xx  = pmax(lwr,pmin(upr,xx))
ProjGrad = pmax( lwr-xx , pmin( upr-xx , -gg ))   
KKTpoint = FALSE
if (max(abs(ProjGrad)) < kktTol)   KKTpoint = TRUE
return(KKTpoint)
}

#------------------------------------------------------------------------------
GetGradient <- function(fname,gname,fc,xc,h,lwr,upr,n,CompareGrad, ...)    {

  gc = NULL
  fcount = 0
  # If a routine for calculating the gradient has been given, use it.
  if (!is.null(gname))  gc = gname(xc, ...)
  if (!is.null(gc) & !CompareGrad) {
    gradlist <- list(gc = gc, fcount = fcount)
    return(gradlist)
  }
  # Otherwise if no analytic gradient, do finite differences.
  xnew = xc  
  gc = rep(0, times = n)
  for (ii in 1:n)    {
    hh = h*max(1,abs(xc[ii]))
    if (xc[ii] < upr[ii]-hh && xc[ii] > lwr[ii]+hh)  {
      # In the centre of the box so use central differences
      xnew[ii] = xc[ii]+hh  
      ffwd = fname(xnew, ...)  
      xnew[ii] = xc[ii]-hh  
      fbak = fname(xnew, ...)  
      gc[ii] = (ffwd - fbak)/(2*hh)
      fcount = fcount + 2
    } else if (xc[ii] < upr[ii]-hh)   {
      # near lower bound so use forward differences
      xnew[ii] = xc[ii]+hh  
      ffwd = fname(xnew, ...)
      gc[ii] = (ffwd - fc)/hh
      fcount = fcount + 1
    } else if (xc[ii] > lwr[ii]+hh)    {
      # near upper bound so use backward differences
      xnew[ii] = xc[ii]-hh  
      fbak = fname(xnew, ...)  
      gc[ii] = (fc - fbak)/hh
      fcount = fcount + 1
    } else {
      gc[ii] = 0
    }
    xnew[ii] = xc[ii]
  }
  if (CompareGrad)   {
    gcc = gname(xc, ...)
    print("finite difference estimate listed first then exact gradient")
    print(gc)
    print(gcc)
  }
  gradlist <- list(gc = gc, fcount = fcount)
  return(gradlist)
}

#------------------------------------------------------------------------------
UpdateB <- function(B,xc,gc,oldxc,oldgc,y,n)    {
  # Update the matrix B using the BFGS update with minor modifications
  s = xc - oldxc
  y = gc - oldgc
  sty = max(1e-6,sum(s*y))
  Bs = as.vector(B%*%s)
  sBs = sum(s*Bs)
  sBsinv = sign(sBs)/max(1e-6,abs(sBs))
  newB = B - outer(sBsinv*Bs,Bs) + outer(y,y/sty)
  return(newB)
}

#------------------------------------------------------------------------------
CheckGradient <- function(Gradient,n)   {
  # Checks the gradient vector for NaN, NA, NULL
  GradOk = TRUE
  for (ii in 1:n)  {
    if (is.null(Gradient[ii]) | !is.finite(Gradient[ii]))  GradOk = FALSE
  }
  return(GradOk)
}

#------------------------------------------------------------------------------
# The QP subproblem solver and its subroutines start here.
#------------------------------------------------------------------------------
CLboxQPsolver <- function(xk,Gradient,Hessian,Delta,upr,lwr,n,QPGradTol)   {
  # Solves the box constrained QP    min_s 0.5 (s^T Hessian s) + s^T Gradient 
  # subject to   lwr <= xk+s <= upr  and   ||Psi^{-1} s|| <= Delta 
  
  # Calculate vector of Coleman-Li scale factors psi and their Jacobian 
  GetPsi.out <- GetPsi(xk,Gradient,upr,lwr,n)
  Psi = GetPsi.out$Psi
  Jacobian = GetPsi.out$J
  
  # Translate the IQP to the new IQP in w with objective
  # Q = 0.5 w^T H w + b^T w    and trust region    ||w|| <= Delta
  b = Psi*Gradient
  H = matrix(0,n,n)
  # H = Psi*Hessian*Psi + diag(Gradient)*Jacobian. Off diagonal entries first.
  for (ii in 1:n)  {
    for (jj in 1:n)  {
      H[ii,jj] = Psi[ii]*Psi[jj]*Hessian[ii,jj]
    }
  }
  # Then do diagonal elements of H.
  for (ii in 1:n)  {
    H[ii,ii] = Psi[ii]*Psi[ii]*Hessian[ii,ii] + Jacobian[ii]*Gradient[ii]
  }
  w = rep(0, times = n)   
  oldp = w    
  gogo = TRUE
  ItnNr = 1
  
  # get the initial gradient
  g = b   
  oldg = g  
  normg = sqrt(sum(g*g))
  # if initial gradient = 0 within tolerance, return initial point (0 step).
  if (normg < QPGradTol)  {  
    CLboxQPsolver.out <- list(xnew = xk, qnew = 0)
    return(CLboxQPsolver.out)    
  }
  
  # Start of the main loop of the QP solver
  while (gogo)   {
    y = g - oldg   
    TRstruck = FALSE   
    # Form the next Polak-Ribiere CG search direction
    p = -g + oldp*(sum(y*g)/(sum(oldg*oldg)))
    # Reset to steepest descent if p not descent direction
    if (sum(p*g) >= 0)  p = -g
    # Calculate the alpha for an exact line search
    denom = sum(p*as.vector(H%*%p))
    if (denom > 0) {
      # positive curvature along search direction p
      alpha = - sum(g*p) / denom
    } else {
      # curvature <= 0 along p, so step to trust region boundary
      TRstruck = TRUE
      wtp = sum(w*p)  
      ptp = sum(p*p)  
      wtw = sum(w*w)
      discriminant = wtp^2 - ptp*(wtw-Delta^2)
      if (discriminant > 0) {
        alpha = (-wtp + sqrt(discriminant) )/ptp
        alphaTR = max(0,alphaTR)
      } else {
        # already at trust region boundary
        alpha = 0  
      }
    }
    
    # Form the proposed step to the new w value: wtest
    wtest = w + alpha*p
    if (sum(wtest*wtest) < Delta*Delta) {
      # step does not reach the trust region boundary TRR
      alphaTR = alpha
    } else {
      # Trust region radius exceeded: shorten step to TRR
      TRstruck = TRUE
      wtp = sum(w*p)  
      ptp = sum(p*p)  
      wtw = sum(w*w)
      discriminant = wtp^2 - ptp*(wtw-Delta^2)
      if (discriminant > 0)  {
        alphaTR = (-wtp + sqrt(discriminant) )/ptp
        alphaTR = max(0,alphaTR)
      } else {
        # already at trust region boundary, so take zero step
        alphaTR = 0  
      }
      wtest = w + alphaTR*p
    }
    # Find maximum zeta such that xk + zeta*s = xk + zeta*Psi*w is feasible
    xs = xk + Psi*wtest 
    ps = Psi*p
    zeta = Inf
    for (ii in 1:n)   {
      if (ps[ii] > 0) {
        zeta = min(zeta,(upr[ii] - xs[ii])/ps[ii])
      } else if (ps[ii] < 0)   {
        zeta = min(zeta,(lwr[ii] - xs[ii])/ps[ii])
      }
    }
    zeta = max(zeta,0)
    # If the trust region boundary has been reached and is feasible, exit
    if (TRstruck | zeta <= alphaTR)    { 
      wtest = w + min(zeta,alphaTR)*p
      sk = Psi*wtest   
      q = sum( sk*(as.vector(Hessian%*%(sk/2)) + Gradient) )
      xnew = pmin(upr,pmax(lwr,xk+sk))
      CLboxQPsolver.out <- list( xnew = xnew,qnew = q)
      return(CLboxQPsolver.out)
    }
    # Take shorter of full step and step along p to edge of feasible region.
    w = w + min(alphaTR,zeta)*p    
    oldp = p 

    ItnNr = ItnNr+1
    oldg = g 
    g = as.vector(H%*%w) + b 
    normg = sqrt(sum(g*g))
    if (normg < QPGradTol)  gogo = FALSE
    if (ItnNr > n)  gogo = FALSE
    sk = Psi*w  
    q = sum( sk*(as.vector(Hessian%*%(sk/2)) + Gradient) )
  }
  # Calculate the step in the original coordinates
  sk = Psi*w 
  q = sum( sk*(as.vector(Hessian%*%(sk/2)) + Gradient) )
  # Project xk + sk onto feasible region and return
  xnew = pmin(upr,pmax(lwr,xk+sk))
  CLboxQPsolver.out <- list( xnew = xnew,qnew = q)
  return(CLboxQPsolver.out)
}

# ----------------------------------------------------------------------
GetPsi <- function(xk,Gradient,upr,lwr,n)   {
  # Calculate the Coleman-Li scale factors Psi and their Jacobian
  # Used in CLboxQPsolver.
  
  xk = pmin(upr,pmax(lwr,xk))
  phi = rep(0, times = n)   
  J = phi
  for (ii in 1:n)   {
    if (Gradient[ii] < 0)   {  
      Theta.out = Theta(upr[ii] - xk[ii])
      xksign = -1
    } else if (Gradient[ii] > 0) {
      Theta.out = Theta(xk[ii] - lwr[ii])
      xksign = 1
    } else {
      Theta.out = Theta( min(upr[ii] - xk[ii],xk[ii] - lwr[ii]) )
      xksign = 1
      if (2*xk[ii] > upr[ii] + lwr[ii])  xksign = -1
    }
    phi[ii] = Theta.out$value
    J[ii] = xksign*Theta.out$deriv
  }
  Psi = sqrt(phi)   
  GetPsi.out <- list(Psi = Psi,J = J)
  return(GetPsi.out)
}

# ----------------------------------------------------------------------
Theta <- function(xi) {
  # Calculate function value and derivative of a single Coleman-Li scale factor
  # Used in CLboxQPsolver.
  
  xi = abs(xi)   
  value = xi/(1+xi) 
  deriv = 1/(1+xi)^2
  Theta.out <- list(value = value,deriv = deriv)
  return(Theta.out)
}

