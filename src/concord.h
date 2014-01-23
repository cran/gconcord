#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

double shrink(double a, double b);
  
int lindx(int r, int c, int p);
  
void concordC(int *nIn, int *pIn, double *S, double *lambda, double *omega, 
	      double *tol, int *maxit, double *iterates);  
  
/* void concordC2(int *nIn, int *pIn, double *S, double *lambda, double *omega, */
/* 	       double *tol, int *maxit, double *iterates); */
void concordC2(int *nIn, int *pIn, double *Y, double *S, double *lambda, double *omega,
	       double *tol, int *maxit, double *iterates);

void symlassoC(int *nIn, int *pIn, double *C, double *lambda, double *omega, 
	       double *tol, int *maxit, double *iterates);
  
void symlassoC2(int *nIn, int *pIn, double *C, double *lambda, double *omega, 
		double *tol, int *maxit, double *iterates);

#endif
