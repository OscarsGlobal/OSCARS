# OSCARS 0.2.1

* `oscars()` gained a `progress` argument (default `TRUE`) that draws a console
  progress bar over the `nfmax` function evaluation budget, showing the
  percentage of that budget used and the elapsed run time. 
  `poscars()` workers run with `progress = FALSE`.
* New import of the `cli` package.

# OSCARS 0.2.0

* Added `poscars()`, a parallel wrapper around `oscars()` to run several
  independent OSCARS-II searches at once (one per processor core) and returns
  the best result. The number of cores is controlled by the new `ncores`
  argument. By default the evaluation budget `nfmax` is divided among the
  workers, cutting wall-clock time by roughly a factor of `ncores`; set
  `divide.budget = FALSE` for a best-of-`ncores` multi-start instead.
* `poscars()` gained a `seed` argument for reproducible parallel runs and a
  `cl` argument for reusing a user-supplied cluster across calls.
* New import of the base `parallel` package.

# OSCARS 0.1.1

* Initial CRAN submission. Initial submission contains 
a single optimization routine, `oscars`. 
