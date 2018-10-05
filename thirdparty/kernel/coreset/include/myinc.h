/*********************************************************
 *
 * myinc.h
 *     - header file for myinc.cpp
 *	 usefur procedures for handling vectors
 *
 * Hai Yu (fishhai@cs.duke.edu), Jun 2003
 *
 \*********************************************************/

#ifndef __MYINC__H
#define __MYINC__H

#define EQUAL(x,y)      (fabs((x)-(y)) < 1E-15)
#define SQR(x)		((x)*(x))

#define pi		3.1415926
#define TIME_INFTY	1E15

/* return squared length of vector */
double v_sqr( double * data, int dim );

/* return length of vector */
double v_len( double * data, int dim );

/* dot product of two vectors */
double v_dot_pdt( double * data1, double * data2, int dim );

/* distance between two points */
double v_dist( double * data1, double * data2, int dim );

/* squared distance between two points */
double v_dist2( double * data1, double * data2, int dim );

/* point-line distance */
double v_pldist( double *p, 			/* the point */
                 double *p1, double *p2, 	/* two points on the line */
                 int dim );			/* dimension */

/* squared point-line distance */
double v_pldist2( double *p, 			/* the point */
                  double *p1, double *p2, 	/* two points on the line */
                  int dim );			/* dimension */

/* subtract two vectors: v1 = v1 - v2 * scalar */
void v_minus( double * data1, double * data2, int dim, double scalar = 1.0 );

/* add two vectors: v1 = v1 + v2 * scalar  */
void v_plus( double * data1, double * data2, int dim, double scalar = 1.0 );

/* normalize a vector */
void v_normalize( double * data, int dim );

#else  /* __MYINC__H */
#error header file myinc.h included twice
#endif /* __MYINC__H */

