/******************************
 * function.h
 * Header file for discrepancy functions
 * 
 * Jeff M Phillips
 * jeffp@cs.duke.edu
 * 8.09.05
 ******************************/

#ifndef FUNCTION__H
#define FUNCTION__H

// for Bernoulli discrepancy function...
#define G 0.01

#include "header.h"

struct discrepancy_function_t {
  double (*f)(double x, double y);
  double (*fx)(double x, double y);
  double (*fy)(double x, double y);
  double (*lambda)(double x, double y);
};
typedef struct discrepancy_function_t d_fxn;

struct approx_discrepancy_function_t {
  int n;              // number hyperplanes
  double *x;          // x-slope of hyperplane
  double *y;          // y-slope of hyperplane
  double *z;          // offset of tangent.
};
typedef struct approx_discrepancy_function_t ad_fxn;


////////// FUNCTIONS /////////////
void get_gridding_fxn (d_fxn *d, ad_fxn *a, double eps, double npts);
void dumb_gridding_fxn (d_fxn *d, ad_fxn *a, double eps, double npts);
void smart_gridding_fxn (d_fxn *d, ad_fxn *a, double eps, double npts);
double eval_approx_fxn (ad_fxn *a, double x, double y);
void set_Poisson_fxn (d_fxn *d);
void set_Gaussian_fxn (d_fxn *d);
void set_Bernoulli_fxn (d_fxn *d);
void set_gamma_fxn (d_fxn *d);


#endif
