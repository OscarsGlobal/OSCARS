# OSCARS-II Global Optimization

An R package that performs black-box global minimization of a general function 
subject to finite bounds on the unknown parameters using a variant of the OSCARS-II 
algorithm (https://doi.org/10.1007/s10898-020-00928-6).  

## Description
Oscars performs black-box minimization (or maximization) of a general function 
subject to bounds on the unknown parameters.   If all bounds
are finite, oscars acts as a global optimization algorithm.  It has been
adapted to handle infinite upper and lower bounds, in which case the method
has the characteristics of a local method for non smooth problems.  Oscars
does not use or assume the existence of derivatives of the objective
function.  It is a low overhead method for cheaply evaluated black-box
functions.

Oscars is a stochastic direct search method which uses only function values
at selected points.   It generates a finite sequence of nested boxes around
a control point, and randomly samples each box once, in turn.   A new set
of nested boxes is formed if the current set is exhausted or a point better
than the control point is found.   In the latter case the better point
replaces the control.   Initially the control point is set to the better
of an internal initial point and a user supplied start point (if given).

From time to time the control is reset alternately to a random point, or
to the best known point.   Each reset marks the end of one cycle and the
start of the next.   All even numbered cycles start with control points
chosen randomly from the feasible region.   All odd numbered cycles (other
than the first) set the control point equal to the best known point.

Oscars either performs a fixed number of function evaluations, or it
halts if progress stalls for a significant period of time.  In both cases it
returns the best known point and the function evaluated at that point.

## Example

To use `oscars`, you must first define the objective function to be minimized.  This 
regular R function 
must take a vector of parameter values as its first argument, and return a scalar.
Additional arguments, beyond the first, can be included. 
The function may return missing (NaN and NA) values. Once the objective function 
is defined, you supply it to `oscars` where you also specify parameter bounds and 
any other (fixed) parameters that are needed in the objective function. 

This is one of the examples in the R documentation files.  It minimizes the 
Branins camel function, which has a minima of *f* = -1.0316 at
(0.0898,0.7127) and (0.0898,-0.7127) with four other local minima.  Here n is the
number of parameters to be minimized with respect to, and lwr and upr are 
the vectors of bounds for the problem.

``` r
# the camel function
camel <- function(par) {
   x = par[1]
   y = par[2]
   f = 4*x^2 - 2.1*x^4 + (1/3)*x^6 + x*y + 4*(y^4-y^2)
   return(f) }
# Minimize in box (-5,5)
out <- oscars(camel, n = 2, lwr = c(-5,-5), upr = c(5,5))
```

Result is
``` r
> out
Minimum value found at -0.0898420652018077, 0.712656434215101
> summary(out)
Minimum value -1.03162845348986 achieved in 1332 evaluations at 
 -0.0898420652018077, 0.712656434215101 
 Evaluations stopped because Optimum function value tolerance reached. 
```

## Installation

Once R is installed, install the current stable release or
the development release with the following commands:

1.  ***Current release***: Install the current release directly from
    CRAN. In the R terminal, issue…

``` r
install.packages("OSCARS")
```

2.  ***Development version***: `oscarsGlobal` is under active development,
    and the ‘development’ version often contains patches between
    official releases. Inspect commit messages for commits following the
    most recent release for a description of the patches. Install the
    development version from
    [GitHub](https://github.com/OscarsGlobal/OSCARS) using:

``` r
if( !require("devtools") ){
  install.packages("devtools")
}
devtools::install_github("OscarsGlobal/OSCARS")
```


### Caveat
Keep the following in mind: Black-box optimization methods for arbitrary 
functions do not and cannot provide certificates of optimality if halted 
after a finite amount of time. 
