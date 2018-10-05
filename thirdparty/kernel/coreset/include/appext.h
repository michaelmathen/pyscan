/*********************************************************\
 *
 * appext.h 
 *     - header file for appext.cpp
 *       computing coresets
 *
 * Hai Yu (fishhai@cs.duke.edu), Jun 2003
 *
 \*********************************************************/

#ifndef __APPEXT__H
#define __APPEXT__H

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>
#include <myinc.h>

typedef double glReal;

class glPoint;
typedef glPoint glVector;

glReal * create_buffer( int size );
glReal * create_buffer( int size, glReal * data );
void release_buffer( glReal * &buf );

class glPoint
{
private:
	int dim;
	glReal * pos;

public:
	void init( int dimension )
	{
		dim = dimension;
		pos = ( glReal * ) malloc( sizeof( glReal ) * dim );
		assert( pos != NULL );
	}

	void init( int dimension,    /* dimension of the point */
	           glReal * data )  /* coordinates of the point */
	{
		init ( dimension );

		for ( int i = 0; i < dim; i++ )
			pos[ i ] = data[ i ];
	}

	void dump()	{ if ( pos != NULL) free( pos ); pos = NULL; }

	glPoint() { dim = 0; pos = NULL; }

	~glPoint() { dump(); }

public:
	int getDim() { return dim; }

	glReal getComponent( int i ) { return pos[ i-1 ]; }  /* 1 <= i <= dim */

	void setComponent( int i, glReal value ) { pos[ i-1 ] = value; }

	glReal getLength()
	{
		glReal l = 0;

		for ( int i = 0; i < dim; i++ )
			l += SQR( pos[i] );

		return sqrt( l );
	}

	void normalize()
	{
		glReal l = getLength();

		if ( l == 0 ) return;
		for ( int i = 0; i < dim; i++ )
			pos[i] /= l;
	}

	void normalize( glReal len )
	{
		glReal l = getLength();

		if ( l == 0 ) return;
		for ( int i = 0; i < dim; i++ )
			pos[i] = pos[i] * len / l;
	}

	glReal dotProduct( glReal * pnt )
	{
		return v_dot_pdt( pnt, pos, dim );
	}

	glReal dotProduct( glVector v )
	{
		return v.dotProduct( pos );
	}
};

class glBBox
{
public:
	int dim;
	glReal * start;
	glReal * vectors;
	
private:
	glReal ratio;
	glReal * unit_len;

public:
	void init(int dimension, glReal * s, glReal * v)
	{
		dim = dimension;
		start = create_buffer( dim, s );
		vectors = create_buffer( dim * dim, v);
		ratio = -1;
		unit_len = NULL;
	}

	void init()
	{
		init( 0, NULL, NULL );
	}

	void dump()
	{
		release_buffer( start );
		release_buffer( vectors );
		release_buffer( unit_len );
		init( 0, NULL, NULL );
	}

	glBBox() { init(); }
	~glBBox() { dump(); }

	/* setup grids in the box. delta defines the size of a single grid */
	/* more precisely, shrinking the bounding box by delta along each dim gives a unit grid */
	void setGrids( glReal delta )
	{
		assert( delta > 0 );
		ratio = delta;
		if ( unit_len == NULL )
			unit_len = create_buffer( dim );
		for ( int i = 0; i < dim; i++ )
			unit_len[ i ] = ratio * v_len( vectors + dim * i, dim );
	}

	void getGridCoordinate( glReal * pnt, int * grid )
	{
		assert ( ratio > 0 );

		glReal * buffer;
		buffer = create_buffer( dim, pnt );
		v_minus( buffer, start, dim );

		for ( int i = 0; i < dim; i++ )
		{
			grid[ i ] = ( int ) floor( v_dot_pdt(buffer, vectors + dim * i, dim )
				 / v_len( vectors + dim * i, dim ) / unit_len[ i ] );

			/* numerical issue */
			if ( grid[ i ] < 0 )
				grid[ i ] = 0;
		}

		release_buffer( buffer );
	}

	void print()
	{
		int i, j;
		printf( "starting vertex: ( " );
		for ( i = 0; i < dim; i++ )
			printf( "%f ", start[i] );
		printf( ")\n" );
		for ( i = 0; i < dim; i++ )
		{
			printf( "vector %d: ( ", i );
			for ( j = 0; j < dim; j++ )
				printf( "%f ", vectors[ i * dim + j ]);
			printf( ")\n" );
		}
	}
};

class glPointSet
{
private:
	glReal * pnt_set;

public:
	int dim;
	int num;
	glReal ** pnt_array;

public:
	void init( int dimension,	/* dimension of the point set */
		       int number)	/* number of points */
	{
		dim = dimension;
		num = number;
		pnt_set = create_buffer( num * dim );
		pnt_array = new glReal*[num];
		for ( number = 0; number < num; number ++ )
			pnt_array[number] = &(pnt_set[number*dim]);
	}

	void init( int dimension,    /* dimension of the point set */
		   int number,       /* total number of points in the set */
		   glReal * data )   /* coordinates of the points */
	{
		init ( dimension, number );
		memcpy( pnt_set, data, sizeof( glReal ) * dim * num );
	}

	void dump()
	{
		release_buffer( pnt_set );
		delete[] pnt_array;
		pnt_array = NULL;
	}

	glPointSet() { dim = num = 0; pnt_set = NULL; }

	~glPointSet() { dump(); }

public:
	void translate( glVector * v )
	{
		int i = 0, j, index = 0;
		while ( i < num )
		{
			j = 0;
			while ( j < dim )
			{
				pnt_set[ index ] += v->getComponent( j+1 );
				index ++;
				j ++ ;
			}
			i ++;
		}
	}

	void rotate( glVector * vs ) /* vs are a number of dim vectors */
	{
		int i, j, k;
		glReal * temp = create_buffer( dim );

		for ( i = 0; i < num; i++ )
		{
			k = i * dim;

			for ( j = 0; j < dim; j++ )
				temp[j] = vs[j].dotProduct( pnt_set + k );

			for ( j = 0; j < dim; j++ )
				pnt_set[ k + j ] = temp[j];
		}

		release_buffer( temp );
	}

	void project( glVector * v, /* projection direction */
		          int pnt )     /* index of the point on the projected hyperplane */
	{
		glReal scalar, * buffer;
		int i, j;

		buffer = create_buffer( dim );

		for ( i = 1; i <= num; i++ )
		{
			getDirection( pnt, i, buffer );
			scalar = - v->dotProduct( buffer ) / SQR( v->getLength() );

			getPoint( i, buffer );
			for ( j = 0; j < dim; j++ )
				buffer[ j ] += scalar * v->getComponent( j+1 );
			setPoint( i, buffer );
		}

		release_buffer( buffer );
	}

	int getNum() { return num; }

	int getDim() { return dim; }

	void getPoint( int i,                /* 1 <= i <= num */
		           glReal * pnt )
	{
		memcpy( pnt, pnt_set + ( i-1 ) * dim, sizeof( glReal ) * dim );
	}

	void setPoint( int i, 	 		/* 1 <= i <= num */
		           glReal * pnt )
	{
		memcpy( pnt_set + ( i-1 ) * dim, pnt, sizeof( glReal ) * dim );
	}

	int getNearestNeighbor( glReal * pt )
	{
		glReal dist, d;
		int nn, i;

		nn = 1;
		dist = v_dist( pnt_set, pt, dim );
		for ( i = 1; i < num;  i++ )
		{
			d = v_dist( pnt_set + i * dim, pt, dim );
			if ( d < dist )
			{
				dist = d;
				nn = i + 1;
			}
		}

		return nn;
	}

	double getExactDiameter();

	/* compute the two extreme points along the specified axis.
	   return the length along the axis and save the indices of the points in i and j */
	glReal getExtremePoints ( int axis,		/* 1 <= axis  <= dim */
		                      int & i, int & j)
	{
		i = 0, j = 0;

		for ( int k = 0; k < num; k++ )
		{
			if ( pnt_set [ k * dim + axis -1 ] < pnt_set [ i * dim + axis -1 ] ) i = k;
			if ( pnt_set [ k * dim + axis -1 ] > pnt_set [ j * dim + axis -1 ] ) j = k;
		}

		i ++, j ++;
		return ( pnt_set [ j * dim + axis -1 ] - pnt_set [ i * dim + axis -1 ] );
	}

	/* compute the two extreme points along the specified direction.
	   return the length along the axis and save the indices of the points in i and j
	   j is farther from the projected hyperplane than i */
	glReal getExtremePoints( glVector * v, int & i, int & j )
	{
		glReal lmin, lmax, l, * buffer;

		buffer = create_buffer( dim );

		i = j = 1;
		getPoint( 1, buffer );
		lmin = lmax = - v->dotProduct( buffer ) / v->getLength();

		for ( int k = 2; k <= num; k++ )
		{
			getPoint( k, buffer );
			l = - v->dotProduct( buffer ) / v->getLength();
			if ( l < lmin ) { i = k; lmin = l;}
			else if ( l > lmax ) { j = k; lmax = l; }
		}

		release_buffer( buffer );
		return ( lmax - lmin );
	}

	glReal getWidth( glVector * v )
	{
		int i, j;
		return getExtremePoints( v, i, j );
	}

	glReal getWidth( glVector *v, int k, glReal *w );

	glReal getDistance ( int i, int j )
	{
		glReal d = 0;
		for ( int k = 0; k < dim; k++)
			d += SQR( pnt_set[ ( i-1 ) * dim + k ] - pnt_set[ ( j-1 ) * dim + k ] );

		return sqrt( d );
	}

	void getDirection( int i, int j, 
		               glReal * v)
	{
		for ( int k = 0; k < dim; k++ )
			v[ k ] = pnt_set[ ( j-1 ) * dim + k ] - pnt_set[ ( i-1 ) * dim + k ];
	}

	void print()
	{
		int i, j;

		for ( i = 0; i < num; i++ )
		{
			printf("point %d: ( ", i+1);
			for ( j = 0; j < dim; j++ )
				printf("%f ", pnt_set[ i * dim + j ]);
			printf(")\n");
		}
	}

	/* estimate the diameter of the point set within factors sqrt(2), sqrt(3)
	   in two and three dimensions respectively and 2 in higher dimensions.
	   save the indices of the two points realizing the estimation in index1 and index2 */
	glReal getApproxDiam( int & index1, int & index2 );

	/* compute approximate minimum-volume bounding box of the point set.
	   algorithm based on the paper
	       Efficiently approximating the minimum-volume bounding box of a
		   point set in three dimensions
	   by Barequet and Har-Peled. */
	glBBox * getApproxMVBB( double ** bdpts = NULL );

	/* compute a coreset by gridding.
	   return a pointer to the glPointSet structure of the coreset */
	glPointSet * getCoreSet( int grids );

	/* compute a coreset by Dudley method.
	   return a pointer to the glPointSet structure of the coreset */
	glPointSet * getCoreSet2( int grids );

	/* compute a coreset by Dudley and ANN */
	glPointSet * getCoreSet3( int grids, double eps, unsigned char * ff = NULL );

	/* best coreset procedure of all. Dudley + ANN + better grids */
	glPointSet * getCoreSet4( int grids, double eps, unsigned char * ff = NULL, double *gridpts = NULL );
				// parameters "ff" and "gridpts" are for getRobustCoreset() only
				// "grids" should be assigned to the total number of grid points (this is 
				//                 unlike previous three procedures)
		

	/* distribute points uniformly on sphere */
	/* NOTE: it is very slow in high dimensions */
	double* uniform_sph( int d, int n );
	
	/* compute robust coreset */
	glPointSet * peelOneLevel( int grids, double eps, double *gridpts = NULL );
	glPointSet * getRobustCoreset( int grids, double eps, int k = 1 );
	bool merge( glPointSet *p );
};

#else  /* __APPEXT__H */
#error header file appext.h included twice
#endif /* __APPEXT__H */

