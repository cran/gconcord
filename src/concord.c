#include <sys/param.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include "concord.h"

inline double shrink(double a, double b) {
  if (b < fabs(a)) {
    if (a > 0) return(a-b);
    else       return(a+b);
  } else {  
    return(0.0);
  }
}

inline int lindx(int r, int c, int p) {

  int rr,cc;

  if (r < c){ rr = r; cc = c; }
  else      { rr = c; cc = r; }

  return( (p*rr)+cc - ((rr*(rr+1))/2) );

}

void concordC(int *nIn, int *pIn, double *S, double *lambda, double *omega,
	      double *tol, int *maxit, double *iterates) {

  int i,j,k,r; 
  int zero = 0;
  int *rptr;
  int p = *pIn;
  int tp = p*(p+1)/2;
  int converged=0;
  int saveiter=0;
  double sSum,s1,s2;
  double maxdiff=1.;

  double *omegaold;
  omegaold = (double *)malloc(tp*sizeof(double));
  if (omegaold == 0){
    Rprintf("Out of Memory!\n");
    return;
  }

  rptr = &zero;
  if (iterates[0]>0) {
    saveiter = 1;
    rptr = &r;
  }
  r = 0;
  while ( (r<*maxit) && !converged ) {
    maxdiff = 0;

    // Update off-diagonals of Omega
    for (i=0;i<(p-1);i++){
      for (j=i+1;j<p;j++){
	s1 = 0;
	s2 = 0;
	for (k=0;k<p;k++){
	  s1 += omega[i*p+k]*S[j*p+k];
	  s2 += omega[k*p+j]*S[i*p+k];
	}
	s1 -= omega[i*p+j]*S[j*p+j];
	s2 -= omega[i*p+j]*S[i*p+i];
	
	sSum = S[i*p+i] + S[j*p+j];

	/* omega[i*p+j] = shrink( -(s1+s2)/sSum, (*lambda)*(*nIn)/sSum ); */
	omega[i*p+j] = shrink( -(s1+s2)/sSum, lambda[lindx(i,j,p)]/sSum );
	omega[j*p+i] = omega[i*p+j];

	maxdiff = fmax2(maxdiff, fabs(omegaold[lindx(i,j,p)]-omega[i*p+j]));
	omegaold[lindx(i,j,p)] = omega[i*p+j];
      }
    }

    // Update diagonals of Omega
    for (i=0;i<p;i++){
      s1 = 0;
      for (k=0;k<p;k++){
	s1 += omega[i*p+k]*S[i*p+k];
      }
      s1 -= omega[i*p+i]*S[i*p+i];

      omega[i*p+i] = (-s1 + sqrt((s1*s1) + (4*S[i*p+i])))/(2*S[i*p+i]);

      maxdiff = fmax2(maxdiff, fabs(omegaold[lindx(i,i,p)]-omega[i*p+i]));
      omegaold[lindx(i,i,p)] = omega[i*p+i];

    }

    // Save iterates
    if (saveiter > 0){
      for (i=0;i<tp;i++){
	iterates[(*rptr)*tp+i] = omegaold[i];
      }
    }

    // check convergence
    if (maxdiff<*tol){
      converged = 1;
    }
    r++;
    
  }
  *maxit = r; // pass it back as number of iterations
  free(omegaold);

}

void printmat(int n, int p, double *X){

  int i,j;
  for (i=0; i<n; i++){
    Rprintf("|");
    for (j=0; j<p; j++){
      Rprintf("% 8.7f ", X[j*n+i]);
    }
    Rprintf("|\n");
  }
  Rprintf("\n");
}

void concordC2(int *nIn, int *pIn, double *Y, double *S, double *lambda, double *omega,
	      double *tol, int *maxit, double *iterates) {

  int i,j,r;
  int zero = 0;
  int one = 1;
  int *rptr;
  int n = *nIn;
  int p = *pIn;
  int np = n*p;
  int tp = p*(p+1)/2;
  int converged=0;
  int saveiter=0;
  double sSum,s1,s2;
  double maxdiff=1.;
  double tmpdiff=1.;

  double *omegaold;
  omegaold = (double *)malloc(tp*sizeof(double));
  if (omegaold == 0){
    Rprintf("Out of Memory!\n");
    return;
  }
  for (i=0; i<tp; i++){ omegaold[i] = 0; }
  for (i=0; i<p;  i++){ omegaold[lindx(i,i,p)] = 1.; }

  rptr = &zero;
  if (iterates[0]>0) {
    saveiter = 1;
    rptr = &r;
  }

  /*********************/
  /* Compute residuals */
  /*********************/
  double *resid;
  resid = (double *)malloc(np*sizeof(double));
  if (omegaold == 0){
    Rprintf("Out of Memory!\n");
    return;
  }
  F77_NAME(dcopy)(&np, Y, &one, resid, &one);

  /* printmat(n,p,resid); */

  r = 0;
  while ( (r<*maxit) && !converged ) {
    maxdiff = 0;

    // Update off-diagonals of Omega
    for (i=0;i<(p-1);i++){
      for (j=i+1;j<p;j++){

	s1 = -omega[i*p+j]*S[j*p+j] + 
	  omega[i*p+i] * F77_NAME(ddot)(&n, &Y[j*n], &one, &resid[i*n], &one)/n;
	s2 = -omega[j*p+i]*S[i*p+i] + 
	  omega[j*p+j] * F77_NAME(ddot)(&n, &Y[i*n], &one, &resid[j*n], &one)/n;
	
	sSum = S[i*p+i] + S[j*p+j];
	
	omega[i*p+j] = shrink( -(s1+s2)/sSum, lambda[lindx(i,j,p)]/sSum );
	omega[j*p+i] = omega[i*p+j];

	tmpdiff = omega[i*p+j] - omegaold[lindx(i,j,p)];

	if (tmpdiff != 0) {
	  // update residuals
	  sSum = tmpdiff/omega[i*p+i];
	  F77_NAME(daxpy)(&n, &sSum, &Y[j*n], &one, &resid[i*n], &one);
	  sSum = tmpdiff/omega[j*p+j];
	  F77_NAME(daxpy)(&n, &sSum, &Y[i*n], &one, &resid[j*n], &one);
	  
	  maxdiff = fmax2(maxdiff, fabs(tmpdiff));
	  omegaold[lindx(i,j,p)] = omega[i*p+j];
	} 
      }
    }
    
    // Update diagonals of Omega
    for (i=0;i<p;i++){

      s1 = -omega[i*p+i]*S[i*p+i] + 
	omega[i*p+i] * F77_NAME(ddot)(&n, &Y[i*n], &one, &resid[i*n], &one)/n;
      omega[i*p+i] = (-s1 + sqrt((s1*s1) + (4*S[i*p+i])))/(2*S[i*p+i]);

      sSum = omegaold[lindx(i,i,p)]/omega[i*p+i];
      F77_NAME(dscal)(&n, &sSum, &resid[i*n], &one);
      sSum = 1 - sSum;
      F77_NAME(daxpy)(&n, &sSum, &Y[i*n], &one, &resid[i*n], &one);

      maxdiff = fmax2(maxdiff, fabs(omegaold[lindx(i,i,p)]-omega[i*p+i]));
      omegaold[lindx(i,i,p)] = omega[i*p+i];
    }

    // Save iterates
    if (saveiter > 0){
      for (i=0;i<tp;i++){
	iterates[(*rptr)*tp+i] = omegaold[i];
      }
    }

    // check convergence
    if (maxdiff<*tol){
      converged = 1;
    }
    r++;
    
  }
  *maxit = r; // pass it back as number of iterations
  free(omegaold);

}

/* void getnz (int *k, int *i, int *j, int *p){ */
/*   k <- k-1; */
/*   i <- floor(((2*p-1) - sqrt((2*p-1)^2-8*k))/2); */
/*   k0 <- i*(2*p-i-1)/2; */
/*   j <- k-k0; */
/* } */

void concordC3(int *nIn, int *pIn, double *S, double *lambda, double *omega,
	      double *tol, int *maxit, double *iterates) {

  int i,j,r; 
  int zero = 0;
  int one = 1;
  int *rptr;
  int p = *pIn;
  int tp = p*(p+1)/2;
  int converged=0;
  int saveiter=0;
  double sSum,s1,s2;
  double maxdiff=1.;

  double *omegaold;
  omegaold = (double *)malloc(tp*sizeof(double));
  int *activeset;
  activeset = (int *)malloc(tp*sizeof(int));
  if (omegaold == 0 || activeset == 0){
    Rprintf("Out of Memory!\n");
    return;
  }
  for (i=0; i<tp; i++){ omegaold[i] = 0.; }
  for (i=0; i<p;  i++){ omegaold[lindx(i,i,p)] = 1.; }

  rptr = &zero;
  if (iterates[0]>0) {
    saveiter = 1;
    rptr = &r;
  }
  r = 0;
  while ( (r<*maxit) && !converged ) {
    maxdiff = 0;

    // Update off-diagonals of Omega
    for (i=0;i<(p-1);i++){
      for (j=i+1;j<p;j++){

	s1 = F77_NAME(ddot)(&p, &omega[i*p], &one, &S[j*p], &one) - omega[i*p+j]*S[j*p+j];
	s2 = F77_NAME(ddot)(&p, &omega[j*p], &one, &S[i*p], &one) - omega[i*p+j]*S[i*p+i];

	/* Rprintf("% 4.3f % 4.3f\n", s1, s2); */

	sSum = S[i*p+i] + S[j*p+j];

	/* omega[i*p+j] = shrink( -(s1+s2)/sSum, (*lambda)*(*nIn)/sSum ); */
	omega[i*p+j] = shrink( -(s1+s2)/sSum, lambda[lindx(i,j,p)]/sSum );
	omega[j*p+i] = omega[i*p+j];

	maxdiff = fmax2(maxdiff, fabs(omegaold[lindx(i,j,p)]-omega[i*p+j]));
	omegaold[lindx(i,j,p)] = omega[i*p+j];
      }
    }

    // Update diagonals of Omega
    for (i=0;i<p;i++){

      s1 = F77_NAME(ddot)(&p, &omega[i*p], &one, &S[i*p], &one) - omega[i*p+i]*S[i*p+i];
      omega[i*p+i] = (-s1 + sqrt((s1*s1) + (4*S[i*p+i])))/(2*S[i*p+i]);

      maxdiff = fmax2(maxdiff, fabs(omegaold[lindx(i,i,p)]-omega[i*p+i]));
      omegaold[lindx(i,i,p)] = omega[i*p+i];

    }

    // Save iterates
    if (saveiter > 0){
      for (i=0;i<tp;i++){
	iterates[(*rptr)*tp+i] = omegaold[i];
      }
    }

    // check convergence
    if (maxdiff<*tol){
      converged = 1;
    }
    r++;
    
  }
  *maxit = r; // pass it back as number of iterations
  free(omegaold);

}
