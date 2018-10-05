/*********************************************************\
 *
 * myinc.cpp
 * 	useful procedures for handling vectors
 *
 * Hai Yu (fishhai@cs.duke.edu), Jun 2003
 *
 \*********************************************************/

#include <math.h>
#include <myinc.h>

double v_sqr( double * data, int dim )
{
	double r = 0;
	for ( int i = 0; i < dim; i++ )
		r += SQR( data[ i ] );

	return r;
}

double v_len( double * data, int dim )
{
	return sqrt( v_sqr( data, dim ) );
}

double v_dot_pdt( double * data1, double * data2, int dim )
{
	double r = 0;
	for ( int i = 0; i < dim; i++ )
		r += data1[ i ] * data2[ i ];

	return r;
}

double v_dist( double * data1, double * data2, int dim )
{
	return sqrt( v_dist2( data1, data2, dim) );
}

double v_dist2( double * data1, double * data2, int dim )
{
	double r = 0;
	for ( int i = 0; i < dim; i++ )
		r += SQR( data1[ i ] - data2[ i ] );

	return r;
}

double v_pldist( double *p, double *p1, double *p2, int dim )
{
	return sqrt( v_pldist2( p, p1, p2, dim ) );
}

double v_pldist2( double *p, double *p1, double *p2, int dim )
{
	double a, b, c;

	a = v_dist2( p, p1, dim );
	b = v_dist2( p1, p2, dim );
	c = v_dist2( p, p2, dim );

	return (4*a*b-SQR(a+b-c))/(4*b);
}

void v_minus( double * data1, double * data2, int dim, double scalar )
{
	for ( int i = 0; i < dim; i ++ )
		data1[ i ] -= data2[ i ] * scalar;
}

void v_plus( double * data1, double * data2, int dim, double scalar )
{
	for ( int i = 0; i < dim; i ++ )
		data1[ i ] += data2[ i ] * scalar;
}

void v_normalize( double * data, int dim )
{
	double l = v_len( data, dim );
	for ( int i = 0; i < dim; i ++ )
		data[i] /= l;
}
