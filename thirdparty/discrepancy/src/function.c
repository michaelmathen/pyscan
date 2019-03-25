/*******************************
 * function.c
 * Deals with functions
 * 
 * Jeff M Phillips
 * jeffp@cs.duke.edu
 * 8.08.05
 ******************************/

#include "function.h"
#include "math.h"

 double dist_sq(double x1, double y1, double x2, double y2) {
  double dx = x1-x2;
  double dy = y1-y2;
  return dx*dx + dy*dy;
}

/**
 * Determines a sparse gridding of the x,y in [0,1] such that < (1+eps) error.
 * 
 * @param d       Discrepancy function to approximate
 * @param a       Approximate discrepancy function
 * @param eps     Acceptable approximation error
 * @param npts    Number of red points 
 */
void get_gridding_fxn (d_fxn *d, ad_fxn *a, double eps, double npts) {
  int i;
  double X,Y;
  int s = 0;
  int c = 32;
  double bd = (double)BASE/npts;
  double l, d2;
  int num_to_add;

  a->x = (double*) malloc (sizeof(double) * c);
  a->y = (double*) malloc (sizeof(double) * c);
  a->z = (double*) malloc (sizeof(double) * c);

  do {
    // calc what to add.
    l = d->lambda(1-bd, bd);
    d2 = sqrt(2*eps/l);
    num_to_add = (1-2*bd)/d2;

    //resize
    if (c<s + num_to_add*2 + 1) { 
      do {c*=2;} while (c<s + num_to_add*2 + 1);
      a->x = (double*) realloc (a->x, sizeof(double) * c);
      a->y = (double*) realloc (a->y, sizeof(double) * c);
      a->z = (double*) realloc (a->z, sizeof(double) * c);
    }

    // calculate hyperplanes
    X = 1-bd;
    for (i=0; i<num_to_add+1; ++i) {
      Y = bd+i*d2;
      a->x[s+i] = d->fx(X,Y);
      a->y[s+i] = d->fy(X,Y);
      a->z[s+i] = d->f(X,Y) - X*a->x[s+i] - Y*a->y[s+i];
    }
    s+=num_to_add;
    Y = bd;
    for (i=1; i<num_to_add+1; ++i) {
      X = 1-bd-i*d2;
      a->x[s+i] = d->fx(X,Y);
      a->y[s+i] = d->fy(X,Y);
      a->z[s+i] = d->f(X,Y) - X*a->x[s+i] - Y*a->y[s+i];
    }
    s+=num_to_add+1;

    bd += d2;
  } while (bd < .5);

  // place one hyperplane at minimum  [CHECK THIS IS always (.5,.5) !!!!]
  X=.5; Y=.5;
  a->x[s] = d->fx(X,Y);
  a->y[s] = d->fy(X,Y);
  a->z[s] = d->f(X,Y) - X*a->x[s] - Y*a->y[s];
  s++;

  // final settings in approx. disc. fxn.
  a->x = (double*) realloc (a->x, sizeof(double) * s);
  a->y = (double*) realloc (a->y, sizeof(double) * s);
  a->z = (double*) realloc (a->z, sizeof(double) * s);  
  a->n = s;
}

/**
 * Determines a uniform gridding of the x,y in [0,1] such that < (1+eps) error.
 * 
 * @param d       Discrepancy function to approximate
 * @param a       Approximate discrepancy function
 * @param eps     Acceptable approximation error
 * @param npts    Number of red points 
 */
void dumb_gridding_fxn (d_fxn *d, ad_fxn *a, double eps, double npts) {
  int i,j,k;
  double X,Y;
  double bd = (double)BASE/npts;
  double l = d->lambda(1-bd, bd);
  double d2 = sqrt(2*eps/l);
  int row_num = (1-2*bd)/d2 + 1;
  double delta = (1-2*bd)/(row_num-1);

  a->x = (double*) malloc (sizeof(double) * (row_num*row_num +1));
  a->y = (double*) malloc (sizeof(double) * (row_num*row_num +1));
  a->z = (double*) malloc (sizeof(double) * (row_num*row_num +1));

  for (i=0; i<row_num; ++i) {
    X = bd+i*delta;
    for (j=0; j<row_num; ++j) {
      Y = bd+j*delta;
      k = i*row_num + j;
      a->x[k] = d->fx(X,Y);
      a->y[k] = d->fy(X,Y);
      a->z[k] = d->f(X,Y) - X*a->x[k] - Y*a->y[k];
    }
  }

  // place one hyperplane at minimum  [CHECK THIS IS always (.5,.5) !!!!]
  X=.5; Y=.5;
  k = row_num*row_num;
  a->x[k] = d->fx(X,Y);
  a->y[k] = d->fy(X,Y);
  a->z[k] = d->f(X,Y) - X*a->x[k] - Y*a->y[k];
  
  a->n = row_num*row_num + 1;
}


/**
 * Determines a sparser gridding of the x,y in [0,1] such that < (1+eps) error.
 * 
 * @param d       Discrepancy function to approximate
 * @param a       Approximate discrepancy function
 * @param eps     Acceptable approximation error
 * @param npts    Number of red points 
 */
void smart_gridding_fxn (d_fxn *d, ad_fxn *a, double eps, double npts) {
  double X,Y,x,y;
  int s = 0;
  int c = 32;
  double bd = (double)BASE/npts;
  double l, d2;

  a->x = (double*) malloc (sizeof(double) * c);
  a->y = (double*) malloc (sizeof(double) * c);
  a->z = (double*) malloc (sizeof(double) * c);

  do {
    // calculate hyperplanes along Y
    X = 1-bd;    Y = bd;
    l = d->lambda(X,Y);
    d2 = sqrt(2*eps/l);
    x = X-d2/2;  y = Y;
    do {      
      l = d->lambda(X,Y);
      d2 = sqrt(2*eps/l);
      y += d2/2;
      if (y > 1-bd)  break;
      
      //resize
      if (c<=s) { 
	c*=2;
	a->x = (double*) realloc (a->x, sizeof(double) * c);
	a->y = (double*) realloc (a->y, sizeof(double) * c);
	a->z = (double*) realloc (a->z, sizeof(double) * c);
      }
      
      a->x[s] = d->fx(x,y);
      a->y[s] = d->fy(x,y);
      a->z[s] = d->f(x,y) - x*a->x[s] - y*a->y[s];
      
      s++;
      Y += d2;   y = Y;
    } while (Y < 1-bd);


    // calculate hyperplanes along X
    X = 1-bd;    Y = bd;
    l = d->lambda(X,Y);
    d2 = sqrt(2*eps/l);
    X -= d2;
    x = X;       y = Y+d2/2;
    while (X > bd) {      
      l = d->lambda(X,Y);
      d2 = sqrt(2*eps/l);
      x -= d2/2;
      if (x < bd)  break;

      //resize
      if (c<=s) { 
	c*=2;
	a->x = (double*) realloc (a->x, sizeof(double) * c);
	a->y = (double*) realloc (a->y, sizeof(double) * c);
	a->z = (double*) realloc (a->z, sizeof(double) * c);
      }
      
      a->x[s] = d->fx(x,y);
      a->y[s] = d->fy(x,y);
      a->z[s] = d->f(x,y) - x*a->x[s] - y*a->y[s];
      
      s++;
      X -= d2;   x = X;
    } 

    X = 1-bd;    Y = bd;
    l = d->lambda(X,Y);  
    d2 = sqrt(2*eps/l);
    bd += d2;
  } while (bd < .5);

  // final settings in approx. disc. fxn.
  a->x = (double*) realloc (a->x, sizeof(double) * (s+1));
  a->y = (double*) realloc (a->y, sizeof(double) * (s+1));
  a->z = (double*) realloc (a->z, sizeof(double) * (s+1));  
  a->n = s+1;

  // place one hyperplane at minimum  [CHECK THIS IS always (.5,.5) !!!!]
  X=.5; Y=.5;
  a->x[s] = d->fx(X,Y);
  a->y[s] = d->fy(X,Y);
  a->z[s] = d->f(X,Y) - X*a->x[s] - Y*a->y[s];
}


/**
 * Evaluates the approximate disc. Function at a point (x,y)
 * 
 * @param a   Approximate discrepancy function - n hyperplanes
 * @param x   x-coord to eval at
 * @param y   y-coord to eval at
 * @return    Approximated discrepancy function
 */
double eval_approx_fxn (ad_fxn *a, double x, double y) {
  int i;
  double m, maxd=0;
  
  for (i=0; i<a->n; ++i) {
    m = x*a->x[i] + y*a->y[i] + a->z[i];
    if (maxd < m) maxd = m;
  }

  return maxd;
}



/**
 * Sets the discrepancy function for points from a Poisson function
 * 
 * @param d   Structure for the discrepancy function
 */
 double Pf(double x, double y) {return x*log10(x/y) + (1-x)*log10((1-x)/(1-y));}
 double Pfx(double x, double y) {return log(x/y) - log((1-x)/(1-y));}
 double Pfy(double x, double y) {return (y-x)/(y*(1-y));}
 double Plambda(double x, double y) {
  double a = 1/(x*(1-x));
  double b = 1/(y*(1-y));
  double c = -x/y - (1-x)/((1-y)*(1-y));
  double apc = a+c;
  return .5*(apc) + .5*sqrt(apc*apc - 4*(a*c-b*b));
}
void set_Poisson_fxn (d_fxn *d) {
  d->f = Pf;
  d->fx = Pfx;
  d->fy = Pfy;
  d->lambda = Plambda;
}



/**
 * Sets the discrepancy function for points from a Gaussian function
 * 
 * @param d   Structure for the discrepancy function
 */
 double Gf(double x, double y) {return (x-y)*(x-y)/(y*(1-y));}
 double Gfx(double x, double y) {return 2*(x-y)/(y*(1-y));}
 double Gfy(double x, double y) 
  {return 2*(y-x)/(y*(1-y)) - (x-y)*(x-y)*(1-2*y)/(y*y*(1-y)*(1-y));}
 double Glambda(double x, double y) {
  double Y = y*(1-y);
  double xmy = x-y;
  double ty = 1-2*y;
  double a = 2/Y;
  double b = -2/Y - 2*xmy*ty/(Y*Y);
  double c = 2/Y + (4*xmy*ty - 2*xmy*xmy)/(Y*Y) + xmy*xmy*ty*ty/(Y*Y*Y);
  double apc = a+c;
  return .5*(apc) + .5*sqrt(apc*apc - 4*(a*c-b*b));
}
void set_Gaussian_fxn (d_fxn *d) {
  d->f = Gf;
  d->fx = Gfx;
  d->fy = Gfy;
  d->lambda = Glambda;
}



/**
 * Sets the discrepancy function for points from a Bernoilli function
 * 
 * @param d   Structure for the discrepancy function
 */
 double Bf(double x, double y) {
  return (x*log(x/y) + (1-x)*log((1-x)/(1-y)) + 
	  (y/G - x)*log(1-G*x/y) + ((1-y)/G - 1+x)*log(1-G*(1-x)/(1-y)));
}
 double Bfx(double x, double y) {
  return log(x/y) - log((1-x)/(1-y)) + log(1-G*(1-x)/(1-y)) - log(1-G*x/y);
}
 double Bfy(double x, double y) {
  return log(1-G*x/y)/G - log(1-G*(1-x)/(1-y))/G;
}
 double Blambda(double x, double y) {
  double a = 1/x + 1/(1-x) + G/((1-y)-G*(1-x)) + G/(y-G*x);
  double b = -1/(y-G*x) - 1/((1-y) - G*(1-x));
  double c = -1/y + 1/(1-y) - G*(1-x)/((1-y)*((1-y)-G*(1-x))) - G*x/(y*(y-G*x));
  double apc = a+c;
  return .5*(apc) + .5*sqrt(apc*apc - 4*(a*c-b*b));
}
void set_Bernoulli_fxn (d_fxn *d) {
  d->f = Bf;
  d->fx = Bfx;
  d->fy = Bfy;
  d->lambda = Blambda;
}


/**
 * Sets the discrepancy function for points from a gamma function
 * 
 * @param d   Structure for the discrepancy function
 */
 double gf(double x, double y) {return x*log(x/y) + (1-x)*log((1-x)/(1-y));}
 double gfx(double x, double y) {return log(x/y) - log((1-x)/(1-y));}
 double gfy(double x, double y) {return (y-x)/(y*(1-y));}
 double glambda(double x, double y) {
  double a = 1/(x*(1-x));
  double b = 1/(y*(1-y));
  double c = -x/y - (1-x)/((1-y)*(1-y));
  double apc = a+c;
  return .5*(apc) + .5*sqrt(apc*apc - 4*(a*c-b*b));
}
void set_gamma_fxn (d_fxn *d) {
  d->f = gf;
  d->fx = gfx;
  d->fy = gfy;
  d->lambda = glambda;
}
