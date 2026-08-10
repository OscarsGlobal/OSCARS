# OSCARS 0.2.0

* Added `poscars()`, a parallel wrapper around `oscars()` that runs several
  independent OSCARS-II searches at once (one per processor core) and returns
  the best result. The number of cores is controlled by the new `ncores`
  argument. By default the evaluation budget `nfmax` is divided among the
  workers, cutting wall-clock time by roughly a factor of `ncores`; set
  `divide.budget = FALSE` for a best-of-`ncores` multi-start instead.
* `poscars()` also gains a `seed` argument for reproducible parallel runs and a
  `cl` argument for reusing a user-supplied cluster across calls.
* New dependency on the base `parallel` package.

# OSCARS 0.1.1

* Initial CRAN submission. Initial submission contains 
a single optimization routine, `oscars`. 
