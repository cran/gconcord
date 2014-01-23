\name{concord}
\alias{concord}
\title{CONvex CORrelation selection methoD}
\usage{
concord(data, lambda, tol = 1e-05, maxit = 100, save.iterates = FALSE,
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
Estimates a sparse inverse covariance matrix from a convex
pseudo-likelihood function with lasso L1 penalty
}
\details{
Implements the CONCORD method by Khare, Oh and Rajaratnam
(2013) http://arxiv.org/abs/1307.5381
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
concord(data,2)
}

