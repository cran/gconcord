\name{symlasso}
\alias{symlasso}
\title{Symmetric Lasso (symlasso)}
\usage{
symlasso(data, lambda, tol = 1e-05, maxit = 100, save.iterates = FALSE,
  ...)
}
\arguments{
  \item{data}{Data matrix with n observations (rows) and p
  variables (columns)}

  \item{lambda}{Penalty parameter}

  \item{tol}{Convergence threshold}

  \item{maxit}{Maximum number of iterations before
  termination}

  \item{save.iterates}{Returns iterates if TRUE}

  \item{...}{ignored}
}
\description{
Estimates a sparse inverse covariance matrix from a
pseudo-likelihood function formulation with L1 penalty on
inverse covariance elements.
}
\details{
Implements the Symmetric Lasso method by Friedman, Hastie
and Tibshirani (2010)
http://statweb.stanford.edu/~tibs/ftp/ggraph.pdf
}
\examples{
library(mvtnorm)

## True omega
omega <- matrix(0,3,3)
omega[1,2] <- omega[2,3] <- 2.1
omega <- t(omega) + omega
diag(omega) <- 3

sigma <- solve(omega)

## Generate data
set.seed(60)
data <- rmvnorm(100, rep(0,3), sigma)

## Solve
symlasso(data,2.1)
}

