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
#' @param n The number of parameters with which \code{fname} is optimized.
#'
#' @param lwr A vector of lower bounds for the parameters of \code{fname}.
#' If a single value is supplied it is used for all lower bounds.  Lower bounds
#' of \code{-Inf} are acceptable.
#'
#' @param upr A vector of upper bounds for the parameters of \code{fname}.
#' If a single value is supplied it is used for all upper bounds.  Upper bounds
#' of \code{Inf} are acceptable and signal a local rather than global search.
#'
#' @param ... Additional parameters supplied to function \code{fname}.
#'
#' @param start An optional start point for the algorithm, passed unchanged to
#' each worker's \code{\link{oscars}} call.  If a single value is provided it is
#' used for all dimensions.  Default is \code{NULL}.
#'
#' @param controls A list of OSCARS control parameters, such as the iteration
#' budget and tolerances. See \code{\link{oscars.control}} for the full list and
#' descriptions.  When \code{divide.budget = TRUE} the \code{nfmax} element is
#' split among the workers internally; the \code{controls} you pass in are not
#' modified.
#'
#' @param ncores The number of processor cores (parallel worker processes) to
#' use.  Must be a positive integer.  If \code{ncores <= 1} the optimization is
#' run serially by calling \code{\link{oscars}} directly (no cluster is
#' created).  If \code{ncores} exceeds the number of cores detected on the
#' machine it is reduced to that number with a warning.  Default is
#' \code{max(1, parallel::detectCores() - 1)}, leaving one core free.
#'
#' @param divide.budget Logical.  If \code{TRUE} (the default) the total
#' evaluation budget \code{controls$nfmax} is divided among the \code{ncores}
#' workers, so each worker performs about \code{nfmax / ncores} evaluations. If
#' \code{FALSE}, every worker is given the full \code{nfmax} budget.
#'
#' @param seed An optional integer used to seed the parallel random number
#' streams, making the whole parallel run reproducible.  Each worker
#' still receives its own independent
#' stream.  Default is \code{NULL}, which draws non-reproducible streams.
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
#' @return An object of class \code{"oscars"}: a list with the best parameters
#' found (\code{par}), the objective value there (\code{value}), the total number
#' of function evaluations over all workers (\code{evaluations}), a convergence
#' code (\code{convergence}), a stopping \code{message}, the \code{controls}
#' used, the number of cores used (\code{ncores}), and the vector of best values
#' from every worker (\code{all.values}).
#'
#'
#' @seealso \code{\link{oscars}} for the underlying algorithm and the meaning of
#' the return fields; \code{\link{oscars.control}} for the control parameters;
#' \code{\link[parallel]{makeCluster}} for supplying your own cluster.
#'
#' @examples
#' \donttest{
#' # Camel function with global minima of f = -1.0316 at
#' # (0.0898, 0.7127) and (0.0898, -0.7127) plus four other local minima.
#' camel <- function(par) {
#'   x <- par[1]
#'   y <- par[2]
#'   4*x^2 - 2.1*x^4 + (1/3)*x^6 + x*y + 4*(y^4 - y^2)
#' }
#' # Run four independent searches in parallel and keep the best.
#' out <- poscars(camel, n = 2, lwr = c(-5, -5), upr = c(5, 5), ncores = 4)
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
#' out <- poscars(camel, 2, -5, 5, ncores = 4, divide.budget = FALSE,
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
                   , ncores = max(1, parallel::detectCores() - 1)
                   , divide.budget = TRUE
                   , seed = NULL
                   , cl = NULL
                   ){

  # -- Validate the number of cores ------------------------------------------
  if (length(ncores) != 1 || is.na(ncores)) {
    stop("'ncores' must be a single positive integer.")
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

  # -- Set up the cluster -----------------------------------------------------
  # A PSOCK cluster works on all platforms (including Windows).  A caller
  # supplied cluster is used as-is and left running; a cluster we create here is
  # stopped again before returning.
  own_cluster <- is.null(cl)
  if (own_cluster) {
    avail <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
    if (!is.na(avail) && ncores > avail) {
      warning(sprintf(
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

  # Trent believes the following is an undocumented feature of
  # parallel::clusterSetRNGStream : If no .Random.seed exists in the master
  # process, parallel::clusterSetRNGStream initiates the same RNG on all
  # workers.  This behavior is not documented, as far as Trent can tell.
  # But, if .Random.seed does exist in the master, parallel::clusterSetRNGStream
  # behaves as documented (i.e., iseed = NULL initiates different RNG on
  # each sequence, while iseed = 123 initiates the same RNG on all working.
  # Our fix: if no seed is provided, generate a single random number
  # which places .Random.seed in the master's work environment.
  if (is.null(seed))  stats::runif(1)

  # Give every worker its own independent, and (if seed is set) reproducible,
  # random number stream.  This is what makes the parallel searches diverge.
  parallel::clusterSetRNGStream(cl, iseed = seed)

  # -- Run the independent searches in parallel -------------------------------
  # oscarsFun is passed explicitly so the workers do not need the OSCARS package
  # to be attached; the serial routine travels with the job.
  oscarsFun <- oscars
  runOne <- function(i, oscarsFun, fname, n, lwr, upr, dots, start, wcontrols) {
    args <- c(list(fname, n, lwr, upr), dots
            , list(start = start, controls = wcontrols))
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
