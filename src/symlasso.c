#define R_BLAS_H

#include <R.h>
#include <sys/param.h>
#include <R_ext/BLAS.h>
#include "concord.h"

void symlassoC(int *nIn, int *pIn, double *C, double *lambda, double *omega, 
	       double *tol, int *maxit, double *iterates){

  int i,j,k,r;
  int zero = 0;
  int *rptr;
  int p = *pIn;
  int tp = p*(p+1)/2;
  int converged=0;
  int saveiter=0;
  double sij,sji;
  double maxdiff=1.;
  double deltaij=0.;

  double *omegaold;
  omegaold = (double *)malloc(tp*sizeof(double));
  if (omegaold == 0){
    Rprintf("Out of Memory!\n");
    return;
  }

  double *Z;
  Z = (double *)malloc(p*p*sizeof(double));
  if (Z == 0){
    Rprintf("Out of Memory!\n");
    return;
  }

  double *q, *sig;
  q = (double *)malloc(p*sizeof(double));
  sig = (double *)malloc(p*sizeof(double));
  if (q == 0 || sig == 0){
    Rprintf("Out of Memory!\n");
    return;
  }

  for (i=0;i<p;i++){
    q[i] = 1;
    sig[i] = 1;
    for (j=0;j<p;j++){
      omega[i*p+j] = 0.;
      Z[i*p+j] = 0.;
    }
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

        sij = C[i*p+j] + Z[i*p+j]*sig[j] - omega[i*p+j]*sig[j]*C[i*p+i];
        sji = C[j*p+i] + Z[j*p+i]*sig[i] - omega[j*p+i]*sig[i]*C[j*p+j];
	
        deltaij = omega[i*p+j];
	omega[i*p+j] = shrink(-(sij+sji),(*lambda))/((C[j*p+j]*sig[i])+(C[i*p+i]*sig[j]));
	omega[j*p+i] = omega[i*p+j];
        deltaij = omega[i*p+j]-deltaij; // delta = omega_new - omega_old;

	if (deltaij != 0){
	  q[j] = q[j] + 2*Z[i*p+j]*deltaij + C[i*p+i]*deltaij*deltaij;
	  q[i] = q[i] + 2*Z[j*p+i]*deltaij + C[j*p+j]*deltaij*deltaij;

	  for (k=0;k<p;k++){
	    Z[k*p+j] = Z[k*p+j]+C[i*p+k]*deltaij;
	    Z[k*p+i] = Z[k*p+i]+C[j*p+k]*deltaij;
	  }
	}

	maxdiff = MAX(maxdiff, fabs(omegaold[lindx(i,j,p)]-omega[i*p+j]));
	omegaold[lindx(i,j,p)] = omega[i*p+j];
      }
    }
    // Update diagonals of Omega
    for (i=0;i<p;i++){
      sig[i] = (-1+sqrt(1+4*q[i]*C[i*p+i]))/(2*q[i]);
      maxdiff = MAX(maxdiff, fabs(omegaold[lindx(i,i,p)]-(sig[i])));
      omegaold[lindx(i,i,p)] = sig[i];
    }

    // Save iterates
    if (saveiter > 0){
      for (i=0;i<tp;i++){
	iterates[(*rptr)*tp+i] = omegaold[i];
      }
    }

    // check convergence
    if (maxdiff<*tol)
      converged = 1;

    r++;

  }  

  // put in the diagonal elements
  for (i=0;i<p;i++){
    omega[i*p+i] = 1/sig[i];
  }

  // pass back number of iterations
  *maxit = r;

  free(omegaold);
  free(Z);
  free(q);
  free(sig);
}

void symlassoC2(int *nIn, int *pIn, double *C, double *lambda, double *omega, 
		  double *tol, int *maxit, double *iterates){

  int i,j,k,r;
  int p = *pIn;
  int tp = p*(p+1)/2;
  int converged=0;
  double sij,sji;
  double maxdiff=1.;
  double deltaij=0.;

  double *Z;
  Z = (double *)malloc(p*p*sizeof(double));
  if (Z == 0){
    Rprintf("Out of Memory!\n");
    return;
  }

  double *q, *sig;
  q = (double *)malloc(p*sizeof(double));
  sig = (double *)malloc(p*sizeof(double));
  if (q == 0 || sig == 0){
    Rprintf("Out of Memory!\n");
    return;
  }


  for (i=0;i<p;i++){
    q[i] = 1;
    sig[i] = 1;
    for (j=0;j<p;j++){
      omega[i*p+j] = 0.;
      Z[i*p+j] = 0.;
    }
  }

  r = 1;
  while ( (r<=(*maxit)-1) && !converged ) {

    // Update off-diagonals of Omega
    for (i=0;i<(p-1);i++){
      for (j=i+1;j<p;j++){

        sij = C[i*p+j] + Z[i*p+j]*sig[j] - omega[i*p+j]*sig[j]*C[i*p+i];
        sji = C[j*p+i] + Z[j*p+i]*sig[i] - omega[j*p+i]*sig[i]*C[j*p+j];
	
        deltaij = omega[i*p+j];
	omega[i*p+j] = shrink(-(sij+sji),(*lambda))/((C[j*p+j]*sig[i])+(C[i*p+i]*sig[j]));
	omega[j*p+i] = omega[i*p+j];
        deltaij = omega[i*p+j]-deltaij; // delta = omega_new - omega_old;

	if (deltaij != 0){
	  q[j] = q[j] + 2*Z[i*p+j]*deltaij + C[i*p+i]*deltaij*deltaij;
	  q[i] = q[i] + 2*Z[j*p+i]*deltaij + C[j*p+j]*deltaij*deltaij;
	  
	  for (k=0;k<p;k++){
	    Z[k*p+j] = Z[k*p+j]+C[i*p+k]*deltaij;
	    Z[k*p+i] = Z[k*p+i]+C[j*p+k]*deltaij;
	  }
	}

	/* Rprintf("%3.2f %3.2f\n",sij,sji); */
	/* return; */

      }
    }

    /* for (i=0;i<p;i++){ */
    /*   for (j=0;j<p;j++){ */
    /* 	Rprintf("%6.5f ",omega[i*p+j]); */
    /*   } */
    /*   Rprintf("\n"); */
    /* } */
    /* Rprintf("\n"); */

    // Update diagonals of Omega
    for (i=0;i<p;i++){
      sig[i] = (-1+sqrt(1+4*q[i]*C[i*p+i]))/(2*q[i]);
    }

    // Save iterates
    maxdiff = 0;
    for (i=0;i<p;i++){

      k = lindx(i,i,p);
      iterates[r*tp+k] = 1/sig[i];

      for (j=i;j<p;j++){
	k = lindx(i,j,p);
	iterates[r*tp+k] = omega[i*p+j];
	maxdiff = MAX(maxdiff,fabs(iterates[(r-1)*tp+k]-iterates[r*tp+k]));
	/* Rprintf("%d : %f %f %e %e\n", r,iterates[(r-1)*tp+k],iterates[r*tp+k],maxdiff,*tol); */
      }
    }

    // check convergence
    if (maxdiff<*tol)
      converged = 1;

    r++;

  }  

  // put in the diagonal elements
  for (i=0;i<p;i++){
    omega[i*p+i] = 1/sig[i];
  }

  // pass back number of iterations
  *maxit = r-1;

}
