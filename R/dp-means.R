## original implementation by Johnny Myles White, distributed via GitHub...
##
dp.means <- function(x, lambda=NULL, max.iterations=100, tolerance=10e-3){ # {{{

  n <- nrow(x)
  k <- 1
  assignments <- rep(1, n)
  mu.x <- mean(x$x)
  mu.y <- mean(x$y)
  converged <- FALSE
  iteration <- 0
  ss.old <- Inf
  ss.new <- Inf

  while (!converged && iteration < max.iterations) { # {{{
    iteration <- iteration + 1
    for (i in 1:n) {
      distances <- rep(NA, k)
      for (j in 1:k) distances[j] <- (x[i,'x']-mu.x[j])^2 + (x[i,'y']-mu.y[j])^2
      if (min(distances) > lambda) {
        k <- k + 1
        assignments[i] <- k
        mu.x[k] <- x[i, 'x']
        mu.y[k] <- x[i, 'y']
      } else {
        assignments[i] <- which(distances == min(distances))
      }
    }
    
    for (j in 1:k) {
      if (length(assignments == j) > 0) {
        mu.x[j] <- mean(x[assignments == j, 'x'])
        mu.y[j] <- mean(x[assignments == j, 'y'])
      }
    }
    
    ss.new <- 0
      
    for (i in 1:n) {
      ss.new <- ss.new + 
                (x[i, 'x'] - mu.x[assignments[i]])^2 + 
                (x[i, 'y'] - mu.y[assignments[i]])^2
    }
    ss.change <- ss.old - ss.new
    if (!is.nan(ss.change) && ss.change < tolerance) converged <- TRUE
  } # }}}
  
  centers <- x.frame(x = mu.x, y = mu.y)
  return(list("centers"=centers, 
              "assignments"=factor(assignments),
              "k"=k, 
              "iterations"=iteration))
} # }}}

## hierarchical version of the above
## FIXME: implement this
##
hdp.means <- function(x, lambda=NULL, max.iterations=100, tolerance=10e-3){ #{{{

  stop('Hierarchical DP-means is not implemented yet')

} # }}}

## choose lambda for above by maximizing SSbg/SSwg across a range of lambdas 
## FIXME: implement this
##
chooseLambda <- function(x, lambdas=seq(1, 100, 1)) { # {{{

} # }}}

## for Dirichlet process mixture models
## FIXME: implement this
##
stickBreakingProcess = function(k, alpha) { # {{{
  betas = rbeta(k, 1, alpha)
  remainder = c(1, cumprod(1 - betas))[1:k]
  return(remainder * betas)
} # }}}

