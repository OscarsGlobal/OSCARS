#' @title Parallel OSCARS-II bound constrained global optimization
#'
#' @description
#' \code{poscars} is a parallel wrapper around \code{\link{oscars}} that runs
#' several independent OSCARS-II searches at the same
#' time, one per processor core, and returning the best result found.
#'
#'
#' @param fname An R function to be minimized (or maximized). It must take a
#' vector of parameter values as its first argument and return a scalar.
#' Additional arguments can be supplied via \code{...}.  Missing (\code{NaN} or
#' \code{NA}) function values are acceptable; they are replaced with \code{Inf}
#' when minimizing (or \code{-Inf} when maximizing).  Because the runs execute
#' in separate R processes, any helper objects or functions that \code{fname}
#' references from the global environment should either be defined inside
#' \code{fname} or passed to it through \code{...} (see \code{Details}).
#'
#' @inheritParams oscars
#'
#' @param ncores The number of processor cores (parallel workers) to
#' use.  Must be a positive integer, or \code{Inf} to use every core detected
#' on the machine.  If \code{ncores <= 1} the optimization is
#' run serially by calling \code{\link{oscars}} directly (no cluster is
#' created).  If \code{ncores} exceeds the number of cores detected on the
#' machine it is reduced to that number and a message is printed.  Default is
#' \code{2} per CRAN policies.
#'
#' @param divide.budget Logical.  If \code{TRUE} (the default) the total
#' evaluation budget \code{controls$nfmax} is divided among the \code{ncores}
#' workers.  In this case, each worker performs a max of
#' \code{ceiling(nfmax / ncores)} evaluations. If
#' \code{FALSE}, every worker is given the full \code{nfmax} budget.
#'
#' @param seed An optional integer used to seed the parallel random number
#' streams. Setting this makes the entire parallel run reproducible.
#' Regardless, each worker receives its own seed to ensure independent
#' OSCARS streams across workers.  Default is \code{NULL}, which draws
#' non-reproducible streams.
#'
#' @param cl An optional pre-existing parallel cluster object created by
#' \code{\link[parallel]{makeCluster}}.  If supplied it is used for the runs and
#' left running on exit (the caller is responsible for stopping it), which avoids
#' cluster start-up cost across repeated calls.  If \code{NULL} (the default) a
#' temporary cluster of \code{ncores} workers is created and stopped
#' automatically.
#'
#' @details
#'
#' OSCARS-II is a stochastic direct search: each run draws random sample points
#' inside a sequence of nested, shrinking boxes and periodically restarts from a
#' random or the best known point (see \code{\link{oscars}} for full details).
#' Each run depends on its own chain of random draws and the
#' inner search loop cannot be split across cores.  However,
#' \emph{independent whole runs} are possible and \code{poscars}
#' launches \code{ncores} independent runs, each seeded with its own
#' reproducible random
#' number stream, and keeps whichever run finds the best objective value.
#'
#' Each worker run stops using the ordinary \code{oscars} stopping rules.
#' OSCARS stopping rules are governed by tolerances \code{fTol} and
#' \code{xTol}
#' (see \code{\link{oscars.control}}).
#'
#' By default, \code{poscars} divides the total evaluation budget
#' \code{nfmax} evenly among
#' the \code{ncores} workers (see \code{divide.budget}).  This keeps the total
#' number of function evaluations roughly the same as a single serial
#' \code{oscars} call, but should execute faster by roughly \code{ncores} times.
#' Set \code{divide.budget = FALSE} to give every worker
#' the full budget, which does more total work but improves the chance of
#' locating the global optimum.
#'
#' The objective function \code{fname}, its \code{...} arguments, the
#' bounds and the controls are sent to each worker.  Any objects
#' that \code{fname} uses from the global workspace are \emph{not}
#' automatically exported to the workers.  To be safe, make \code{fname}
#' self-contained or pass everything it needs through \code{...}.
#'
#'
#' @return A list of class \code{"oscars"} containing results of the optimization.
#' This list consists of the following components:
#'   \itemize{
#'      \item \code{par}: vector containing the best known parameters
#'      \item \code{value}: The minimized (or maximized) function value
#'      \item \code{evaluations}: The number of function evaluations used
#'      \item \code{cycles}: The number of cycles used
#'      \item \code{convergence}: 0 if the function and parameter tolerances have been reached;
#'       1 if tolerances have not been reached but function evaluation budget has been exhausted;
#'       2 if bounds are inconsistent.
#'      \item \code{message}: A text string explaining the value in \code{convergence}.
#'      \item \code{controls}: The values of the controls provided to oscars.
#'      \item \code{ncores}: The number of cores used.
#'      \item \code{all.values}: A vector containing the best values returned
#'      by every worker.
#'   }
#'
#' @seealso \code{\link{oscars}} for the underlying algorithm and the meaning of
#' the return fields; \code{\link{oscars.control}} for the control parameters;
#' \code{\link[parallel]{makeCluster}} for supplying your own cluster.
#'
#' @examples
#' \donttest{
#'
#' # Per CRAN policies, these examples run on only 2 cores.
#' # Set ncores = Inf to utilize all cores.
#' # Setting ncores = parallel::detectCores()-1 is a good choice.
#'
#' # Camel function with global minima of f = -1.0316 at
#' # (0.0898, 0.7127) and (0.0898, -0.7127) plus four other local minima.
#' camel <- function(par) {
#'   x <- par[1]
#'   y <- par[2]
#'   4*x^2 - 2.1*x^4 + (1/3)*x^6 + x*y + 4*(y^4 - y^2)
#' }
#' # Run four independent searches in parallel and keep the best.
#' out <- poscars(camel, n = 2, lwr = c(-5, -5), upr = c(5, 5), ncores = 2)
#' out
#'
#' # Reproducible parallel run via the seed argument.
#' out1 <- poscars(camel, 2, -5, 5, ncores = 2, seed = 42)
#' out2 <- poscars(camel, 2, -5, 5, ncores = 2, seed = 42)
#' identical(out1$value, out2$value)
#'
#' # Passing extra arguments to the objective function.
#' # Rosenbrock's "banana" function, global minimum 0 at (a, a^2).
#' rosenbrock <- function(par, a = 1, b = 100) {
#'   (a - par[1])^2 + b*(par[2] - par[1]^2)^2
#' }
#' out <- poscars(rosenbrock, 2, -3, 3, a = 0.5, ncores = 2)
#'
#' # Best-of-ncores multi-start: give every worker the full budget instead
#' # of dividing it, trading more total work for a more thorough search.
#' out <- poscars(camel, 2, -5, 5, ncores = 2, divide.budget = FALSE,
#'                controls = oscars.control(nfmax = 20000, infol = 0))
#'
#' # Reuse a single cluster across several calls to avoid start-up cost.
#' cl <- parallel::makeCluster(2)
#' o1 <- poscars(camel, 2, -5, 5, cl = cl)
#' o2 <- poscars(rosenbrock, 2, -3, 3, a = 0.5, cl = cl)
#' parallel::stopCluster(cl)
#' }
#'
#' @export
poscars <- function(fname
                   , n
                   , lwr
                   , upr
                   , ...
                   , start = NULL
                   , controls = oscars.control()
                   , ncores = 2
                   , divide.budget = TRUE
                   , seed = NULL
                   , cl = NULL
                   ){

  # -- Validate the number of cores ------------------------------------------
  if (length(ncores) != 1 || is.na(ncores)) {
    stop("'ncores' must be a single positive integer.")
  }
  # ncores = Inf asks for every core on the machine.  Resolve it here because
  # as.integer(Inf) is NA, which would derail every test made below.
  if (is.infinite(ncores)) {
    avail <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
    ncores <- if (is.na(avail))  1L  else  max(1L, as.integer(avail))
  }
  ncores <- as.integer(ncores)

  # Capture the extra arguments destined for fname once, up front.
  dots <- list(...)

  # -- Serial fall-back -------------------------------------------------------
  # With a single core (and no user supplied cluster) there is nothing to
  # parallelize, so call oscars directly and avoid all cluster overhead.
  if (ncores <= 1 && is.null(cl)) {
    serial <- do.call(oscars
                    , c(list(fname, n, lwr, upr), dots
                      , list(start = start, controls = controls)))
    serial$ncores <- 1L
    serial$all.values <- serial$value
    return(serial)
  }

  # -- Decide the per-worker evaluation budget --------------------------------
  # Each worker keeps the ordinary oscars stopping rules; we only adjust the
  # maximum-evaluations budget when dividing it among the workers.
  wcontrols <- controls
  if (isTRUE(divide.budget)) {
    wcontrols$nfmax <- max(1L, as.integer(ceiling(controls$nfmax / ncores)))
  }
  # Workers run silently; printing from several processes at once is unreadable.
  wcontrols$infol <- 0
  # A progress bar per worker would be just as unreadable, so the workers
  # get progress = FALSE below.

  # -- Set up the cluster -----------------------------------------------------
  # A PSOCK cluster works on all platforms (including Windows).  A caller
  # supplied cluster is used as-is and left running; a cluster we create here is
  # stopped again before returning.
  own_cluster <- is.null(cl)
  if (own_cluster) {
    avail <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
    if (!is.na(avail) && ncores > avail) {
      message(sprintf(
        "Requested ncores = %d exceeds the %d cores detected; using %d.",
        ncores, avail, avail))
      ncores <- as.integer(avail)
    }
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    # Match the reported core count to the cluster the caller supplied.
    ncores <- length(cl)
  }

  # Give every worker its own independent, and (if seed is set) reproducible,
  # random number stream.  This is what makes the parallel searches diverge.
  parallel::clusterSetRNGStream(cl, iseed = seed)

  # -- Run the independent searches in parallel -------------------------------
  # oscarsFun is passed explicitly so the workers do not need the OSCARS package
  # to be attached; the serial routine travels with the job.
  oscarsFun <- oscars
  runOne <- function(i, oscarsFun, fname, n, lwr, upr, dots, start, wcontrols) {
    args <- c(list(fname, n, lwr, upr), dots
            , list(start = start, controls = wcontrols, progress = FALSE))
    do.call(oscarsFun, args)
  }

  results <- parallel::parLapply(cl, seq_len(ncores), runOne
                               , oscarsFun = oscarsFun
                               , fname = fname, n = n, lwr = lwr, upr = upr
                               , dots = dots, start = start
                               , wcontrols = wcontrols)

  # -- Combine the results ----------------------------------------------------
  # Each element is an "oscars" object.  value already carries the correct sign
  # (oscars flips maximization back before returning), so "best" is the largest
  # value when maximizing and the smallest when minimizing.
  vals <- vapply(results, function(r) {
    v <- r$value
    if (is.null(v) || length(v) != 1) NA_real_ else as.numeric(v)
  }, numeric(1))

  DoMax <- isTRUE(controls$DoMax)
  if (all(is.na(vals))) {
    # Every run failed (e.g. inconsistent bounds); return the first result.
    best_i <- 1L
  } else if (DoMax) {
    best_i <- which.max(vals)
  } else {
    best_i <- which.min(vals)
  }
  best <- results[[best_i]]

  # Total function evaluations is the sum across all workers.
  total_evals <- sum(vapply(results, function(r) {
    e <- r$evaluations
    if (is.null(e)) 0L else as.integer(e)
  }, integer(1)))

  # -- Assemble the returned object, matching oscars' structure ---------------
  solution <- list(par = best$par
                 , value = best$value
                 , evaluations = total_evals
                 , convergence = best$convergence
                 , message = best$message
                 , controls = controls
                 , ncores = ncores
                 , all.values = vals)
  class(solution) <- "oscars"

  # Optional summary line, printed once from the master process.
  if (!is.null(controls$infol) && controls$infol > 0) {
    if (DoMax) {
      cat(sprintf("\n poscars: best Max f = %12.6g over %d parallel runs",
                  solution$value, ncores))
    } else {
      cat(sprintf("\n poscars: best Min f = %12.6g over %d parallel runs",
                  solution$value, ncores))
    }
    cat(sprintf(" using %d total function evaluations.\n\n", total_evals))
  }

  return(solution)
}   # end of poscars
