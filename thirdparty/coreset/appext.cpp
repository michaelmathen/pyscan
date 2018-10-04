/*********************************************************\
 *
 * appext.cpp 
 *     - computing coresets
 *
 * Hai Yu (fishhai@cs.duke.edu), Jun 2003
 *
\*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include <sys/time.h>
#include <ANN/ANN.h>

#include <appext.h>

#define PI 3.1415926535

glReal * create_buffer( int size )
{
	if ( size <= 0 )
		return NULL;

	glReal * buffer;
	buffer = ( glReal * ) malloc( sizeof( glReal ) * size );
	assert( buffer != NULL );
	memset( buffer, 0, sizeof( glReal ) * size );

	return buffer;
}

glReal * create_buffer( int size,
		    glReal * data )
{
	if ( size <= 0 || data == NULL )
		return NULL;

	glReal * buffer;
	buffer = create_buffer( size );
	memcpy( buffer, data, sizeof( glReal ) * size );

	return buffer;
}

void release_buffer( glReal * &buf )
{
	if ( buf != NULL ) 
		free( buf );
	buf = NULL;
}

glReal glPointSet::getApproxDiam( int & index1, int & index2 )
{
	int ind1, ind2, i;
	glReal l = -1.0, ltemp;
	
	for ( i = 1; i <= dim; i++)
	{
		ltemp = getExtremePoints( i, ind1, ind2 );
		if ( ltemp > l)
		{
			index1 = ind1;
			index2 = ind2;
			l = ltemp;
		}
	}

	l = -1 ;
	if ( dim >= 4 )
	{
		for ( i = 1; i <= num; i++ )
		{
			ltemp = getDistance( index1, i );
			if ( ltemp > l )
			{
				index2 = i;
				l = ltemp;
			}
		}
	}

	return l;
}

glBBox * glPointSet::getApproxMVBB( double ** bdpts )
{
	glReal * temp, * vectors, * start_pnt, * buffer;
	glReal l, ltemp;
	int i, j, index, temp1, temp2;
	glVector v;
	glBBox * bbox;

	buffer = create_buffer( dim );

	/* save original data */
	temp = create_buffer( num * dim );
	memcpy( temp, pnt_set, sizeof( glReal ) * num * dim );
	
	vectors = create_buffer( dim * dim );
	start_pnt = create_buffer( dim );

	for ( i = dim; i >= 1 ; i-- )
	{
		l = 0; index = 1;
		for ( j = 2; j <= num; j++ )
		{
			ltemp = getDistance(1, j);
			if ( ltemp > l ) { l = ltemp; index = j; }
		}

		getDirection( index, 1, buffer );
		v.init( dim, buffer );
		v.normalize();
		l = getExtremePoints( & v, temp1, temp2 );
		if ( bdpts != NULL )
		{
			memcpy( bdpts[2*i-1], temp+temp1*dim-dim, sizeof(double)*dim );
			memcpy( bdpts[2*i-2], temp+temp2*dim-dim, sizeof(double)*dim );
		}
		for ( j = 0; j < dim; j++ )
			vectors[ (i-1) * dim + j ] = - l * v.getComponent( j+1 );
		project( & v , temp1 );
		if ( i == 1 ) getPoint( temp1, start_pnt );
		v.dump();
	}
	
	/* restore original data */
	memcpy( pnt_set, temp, sizeof( glReal ) * num * dim );
	release_buffer( temp );

	bbox = ( glBBox * ) malloc( sizeof( glBBox ) );
	assert( bbox != NULL );
	bbox->init(dim, start_pnt, vectors);

	release_buffer( start_pnt );
	release_buffer( vectors );
	release_buffer( buffer );
	return bbox;
}


/* constructing core-sets using grid method */
glPointSet * glPointSet::getCoreSet( int grids )
{
	glBBox * bbox;
	bbox = getApproxMVBB();

	/* compute alpha */

	/* setup grids on the bounding box */
	glReal delta;
	delta = 1.0 / grids;
	bbox->setGrids( delta );

	/* compute the core-set */
	int size, index, i, j, k, * flr, * cl, * gd;
	glReal * buffer;
	unsigned char * flag;

	grids ++;  // boundary effects
	size = 1;
	for ( i = 0; i < dim-1; i ++) size *= grids;
	flr = ( int * ) malloc( sizeof( int ) * size );
	assert( flr != NULL );
	memset( flr, 0, sizeof( int ) * size );
	cl = ( int * ) malloc( sizeof( int ) * size );
	assert( cl != NULL );
	memset( cl, 0, sizeof( int ) * size );
	gd = ( int * ) malloc( sizeof( int ) * dim );
	assert ( gd != NULL );

	buffer = create_buffer( dim );

	flag = ( unsigned char * ) malloc( sizeof ( unsigned char ) * num );
	assert( flag != NULL );
	memset( flag, 0 , sizeof ( unsigned char ) * num );

	for ( i = 1; i <= num; i++ )
	{
		getPoint( i, buffer );
		bbox->getGridCoordinate( buffer, gd );
		
		/* compute the index from gd */
		index = 0;
		for ( j =1, k = 0; k < dim - 1; k++, j *= grids )
			index += j * gd[ k ];

		k = gd[ dim-1 ] + 1;
		if ( flr[ index ] == 0 )
		{
			flr[ index ] = cl[ index ] = i;
			flag[ i-1 ] = 1;
		}
		else
		{
			getPoint( flr[ index ], buffer );
			bbox->getGridCoordinate( buffer, gd );
			if ( k > gd[ dim-1 ] + 1 )
			{
				flag[ i-1 ] = 1;
				if ( cl[ index ] != flr[ index ] )
					flag[ flr[ index ] - 1 ] = 0;
				flr[ index ] = i;
			}

			getPoint( cl[ index ], buffer );
			bbox->getGridCoordinate( buffer, gd );
			if ( k < gd[ dim-1 ] + 1)
			{
				flag[ i-1 ] = 1;
				if ( cl[ index ] != flr[ index ] )
					flag[ cl[ index ] - 1 ] = 0;
				cl[ index ] = i;
			}
		}
	}

	/* create glPointSet for the computed coreset */
	for ( k = 0, i = 0; i < num; i++ ) k += flag[ i ];  // k: total number of pts in core-set

	glPointSet * p;
	p = ( glPointSet * ) malloc( sizeof( glPointSet ) );
	assert( p != NULL );
	p->init( dim, k );
	for ( k = 1, i = 0; i < num; i++ )
		if ( flag[ i ] )
		{
			p->setPoint( k, pnt_set + i * dim );
			k ++;
		}

	free( flr );
	free( cl );
	free( gd );
	free( buffer );
	free( flag );
	bbox->dump();

	return p;
}

/* constructing core-sets using dudley method */
glPointSet * glPointSet::getCoreSet2( int grids )
{
	int i, j, k;

	glPointSet p;
	p.init( dim, num, pnt_set );

	unsigned char * flag;
	flag = ( unsigned char * ) malloc( sizeof ( unsigned char ) * p.getNum() );
	assert( flag != NULL );
	memset( flag, 0 , sizeof ( unsigned char ) * num );

	glBBox * bbox;
	bbox = getApproxMVBB();

	/* transform the point set into unit square */
	glReal * buffer = create_buffer( dim );
	v_minus( buffer, bbox->start, dim );
	glPoint vec;
	vec.init( dim, buffer );
	p.translate( &vec );

	glPoint * vecs;
	vecs = new glPoint[dim];
	for ( i = 0; i < dim; i++ )
	{
		vecs[i].init( dim, bbox->vectors + i * dim );
		vecs[i].normalize( 2 / vecs[i].getLength() );
	}
	p.rotate( vecs );

	for ( i = 1; i <= dim; i++ ) vec.setComponent( i, -1 );
	p.translate( &vec );

	/* create grid points on the sphere and find in p the nearest neighbor of each grid point */
	int * index = new int[ dim -1 ];

	for ( i = 0; i < dim - 1; i ++ )
		index[ i ] = 0;
	int n = 1;
	for ( i = 0; i < dim - 1; i ++ )
		n *= grids;
	for ( i = 0; i < n; i ++ )
	{
		vec.setComponent( dim, 1 );
		for ( j = 0; j < dim -1; j ++ )
			vec.setComponent( j+1, 1/tan( pi/(2 * grids) + pi * index[j] / grids) );
		vec.normalize( sqrt( dim ) + 2 );

		for ( j = 0; j < dim; j ++ )
			buffer[ j ] = vec.getComponent( j+1 );
		flag[ p.getNearestNeighbor( buffer ) - 1 ] = 1;

		for ( j = 0; j < dim; j ++ )
			buffer[ j ] = - buffer[ j ];
		flag[ p.getNearestNeighbor( buffer ) - 1 ] = 1;

		index[ 0 ] ++;
		j = 0;
		while ( (j < dim - 1) && (index[ j ] == grids) )
		{
			index[ j ] = 0;
			j ++;
			if ( j < dim -1 )
				index[ j ] ++;
		}
	}

	/* create glPointSet for the computed coreset */
	for ( k = 0, i = 0; i < p.getNum(); i++ )
		k += flag[ i ];  // k: total number of pts in core-set

	glPointSet * pr;
	pr = ( glPointSet * ) malloc( sizeof( glPointSet ) );
	assert( pr != NULL );
	pr->init( dim, k );
	for ( k = 1, i = 0; i < num; i++ )
		if ( flag[ i ] )
		{
			getPoint( i + 1, buffer );
			pr->setPoint( k, buffer );
			k ++;
		}

	/* wrap up */
	delete[] vecs;
	delete[] index;
	release_buffer( buffer );
	free( flag );
	bbox->dump();
	p.dump();

	return pr;
}

/* constructing core-sets using dudley method and ANN structure */
glPointSet * glPointSet::getCoreSet3( int grids, double eps, unsigned char * ff )
{
	int i, j, k;
	struct timeval start, end;

	glPointSet p;
	p.init( dim, num, pnt_set );

	unsigned char * flag;
	if ( ff == NULL )
	{
		flag = ( unsigned char * ) malloc( sizeof ( unsigned char ) * p.getNum() );
		assert( flag != NULL );
		memset( flag, 0 , sizeof ( unsigned char ) * num );
	}
	else
		flag = ff;

	gettimeofday( &start, (struct timezone *)0 );

	glBBox * bbox;
	bbox = getApproxMVBB();

	/* transform the point set into unit square */
	glReal * buffer = create_buffer( dim );
	v_minus( buffer, bbox->start, dim );
	glPoint vec;
	vec.init( dim, buffer );
	p.translate( &vec );

	glPoint * vecs;
	vecs = new glPoint[dim];
	for ( i = 0; i < dim; i++ )
	{
		vecs[i].init( dim, bbox->vectors + i * dim );
		vecs[i].normalize( 2 / vecs[i].getLength() );
	}
	p.rotate( vecs );

	for ( i = 1; i <= dim; i++ ) vec.setComponent( i, -1 );
	p.translate( &vec );

	/* create grid points on the sphere and find in p the nearest neighbor of each grid point */
	// first init ANN data structure

	ANNkd_tree	*the_tree;
	ANNidxArray	nn_idx;
	ANNdistArray	dists;

	nn_idx = new ANNidx[1];
	dists = new ANNdist[1];

	the_tree = new ANNkd_tree( p.pnt_array, p.num, p.dim );

	gettimeofday( &end, (struct timezone *)0 );

	gettimeofday( &start, (struct timezone *)0 );

	int * index = new int[ dim -1 ];
	for ( i = 0; i < dim - 1; i ++ )
		index[ i ] = 0;
	int n = 1;
	for ( i = 0; i < dim - 1; i ++ )
		n *= grids;
	for ( i = 0; i < n; i ++ )
	{
		vec.setComponent( dim, 1 );
		for ( j = 0; j < dim -1; j ++ )
			vec.setComponent( j+1, tan( pi/(2 * grids) + pi * index[j] / grids) );
		vec.normalize( sqrt( dim ) + 2 );

		for ( j = 0; j < dim; j ++ )
			buffer[ j ] = vec.getComponent( j+1 );
		the_tree->annkSearch( buffer, 1, nn_idx, dists, eps );
		flag[ nn_idx[0] ] = 1;

		for ( j = 0; j < dim; j ++ )
			buffer[ j ] = - buffer[ j ];
		the_tree->annkSearch( buffer, 1, nn_idx, dists, eps );
		flag[ nn_idx[0] ] = 1;

		index[ 0 ] ++;
		j = 0;
		while ( (j < dim - 1) && (index[ j ] == grids) )
		{
			index[ j ] = 0;
			j ++;
			if ( j < dim -1 )
				index[ j ] ++;
		}
	}
	gettimeofday( &end, (struct timezone *)0 );

	/* create glPointSet for the computed coreset */
	for ( k = 0, i = 0; i < p.getNum(); i++ )
		k += flag[ i ];  // k: total number of pts in core-set

	glPointSet * pr;
	pr = ( glPointSet * ) malloc( sizeof( glPointSet ) );
	assert( pr != NULL );
	pr->init( dim, k );
	for ( k = 1, i = 0; i < num; i++ )
		if ( flag[ i ] )
		{
			getPoint( i + 1, buffer );
			pr->setPoint( k, buffer );
			k ++;
		}

	/* wrap up */
	delete the_tree;
	delete[] nn_idx;
	delete[] dists;
	delete[] vecs;
	delete[] index;
	release_buffer( buffer );
	if ( ff == NULL ) free( flag );
	bbox->dump();
	p.dump();

	return pr;
}

/* constructing core-sets using dudley method, ANN structure, and enhanced grid construction */
glPointSet * glPointSet::getCoreSet4( int grids, double eps, unsigned char * ff, double *gridpts )
{
	int i, j, k;
	struct timeval start, end;

	glPointSet p;
	p.init( dim, num, pnt_set );

	/* initializing */
	unsigned char * flag;
	if ( ff == NULL )
	{
		flag = ( unsigned char * ) malloc( sizeof ( unsigned char ) * p.getNum() );
		assert( flag != NULL );
		memset( flag, 0 , sizeof ( unsigned char ) * num );
	}
	else
		flag = ff;

	glBBox * bbox;
	/* compute approximate bounding box */
	bbox = getApproxMVBB();

	/* transform the point set into unit square */
	glReal * buffer = create_buffer( dim );
	v_minus( buffer, bbox->start, dim );
	glPoint vec;
	vec.init( dim, buffer );
	p.translate( &vec );

	glPoint * vecs;
	vecs = new glPoint[dim];
	for ( i = 0; i < dim; i++ )
	{
		vecs[i].init( dim, bbox->vectors + i * dim );
		vecs[i].normalize( 2 / vecs[i].getLength() );
	}
	p.rotate( vecs );

	for ( i = 1; i <= dim; i++ ) vec.setComponent( i, -1 );
	p.translate( &vec );

	/* create grid points on the sphere and find in p the nearest neighbor of each grid point */
	// first init ANN data structure
	ANNkd_tree	*the_tree;
	ANNidxArray	nn_idx;
	ANNdistArray	dists;

	nn_idx = new ANNidx[1];
	dists = new ANNdist[1];

	the_tree = new ANNkd_tree( p.pnt_array, p.num, p.dim );

	// create grid pts
		
	double * grid_pts;
	if ( gridpts != NULL )
		grid_pts = gridpts;
	else
	{
		grid_pts = uniform_sph( dim, grids );
		double r = sqrt(dim) + 2;
		for ( i = 0; i < dim*grids; i ++ )
			grid_pts[i] *= r;
	}

	for ( i = 0; i < grids; i ++ )
	{
		the_tree->annkSearch( &(grid_pts[i*dim]), 1, nn_idx, dists, eps );
		flag[ nn_idx[0] ] = 1;
	}
	if ( gridpts == NULL ) delete[] grid_pts;

	/* create glPointSet for the computed coreset */
	for ( k = 0, i = 0; i < p.getNum(); i++ )
		k += flag[ i ];  // k: total number of pts in core-set

	glPointSet * pr;
	pr = ( glPointSet * ) malloc( sizeof( glPointSet ) );
	assert( pr != NULL );
	pr->init( dim, k );
	for ( k = 1, i = 0; i < num; i++ )
		if ( flag[ i ] )
		{
			getPoint( i + 1, buffer );
			pr->setPoint( k, buffer );
			k ++;
		}

	/* wrap up */
	delete the_tree;
	delete[] nn_idx;
	delete[] dists;
	delete[] vecs;
	release_buffer( buffer );
	if ( ff == NULL ) free( flag );
	bbox->dump();
	p.dump();

	return pr;
}

double glPointSet::getExactDiameter()
{
	double diam, temp;
	int i, j;

	diam = 0;
	for ( i = 0; i < num-1; i ++ )
	{
		if ( i % 100 == 0 ) 
			printf("\r - %.3f %%", i*100.0/num);
		for ( j = i + 1; j < num; j ++ )
		{
			temp = v_dist2( pnt_array[i], pnt_array[j], dim );
			if ( temp > diam )
				diam = temp;
		}
	}
	printf("\n");

	return sqrt(diam);
}

/* distributing n points on (dim-1) unit sphere as uniform as possible */
/* NOTE: this procedure produces very good point distributions, but is 
   extremely slow in high dimensions. I'll probably improve this procedure
   in a future version */
double* glPointSet::uniform_sph( int d, int n )
{
	int i,j,k,m,z;

	// initialize
	double* gridpts_buf = new double[ d*n ];
	double** grid_pts = new double*[ n ];
	for ( i = 0; i < n; i ++ )
		grid_pts[i] = &(gridpts_buf[i*d]);

	// start with uniform distribution
	for ( i = 0; i < n*d; i++ )
		gridpts_buf[i] = 1.0;
	for ( i = 0; i < n; i ++ )
	{
		while ( v_len(grid_pts[i],d) > 1 )
			for ( j = 0; j < d; j++ )
				grid_pts[i][j] = (rand()%100003-50001)/50001.0;

		v_normalize(grid_pts[i], d);
	}

	// gradient descent
	double rate = 0.1;
	double *x = new double[d-1];
	double *dx = new double[d-1];
	double *tmp = new double[d];
	double f, dist;
	for ( i = 0; i < 20; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			//compute x_1,...,x_{d-1}
			f = 1;
			for ( k = 0; k < d-1; k++ )
			{
				x[k] = asin( grid_pts[j][k]/f );
				f *= cos( x[k] );
			}
			if ( grid_pts[j][d-1] < 0 )
				x[d-2] = PI - x[d-2];

			for ( k = 0; k < d-1; k++ )
				dx[k] = 0;
			for ( k = 0; k < n; k++ )
			{
				if ( k == j ) continue;

				dist = v_dist( grid_pts[j], grid_pts[k], d );
				for ( m = 0; m < d-1; m++ )
				{
					f = 1;
					for ( z = 0; z < d-1; z++ )
					{
						if ( z < m )
						{
							tmp[z] = 0;
							f *= cos(x[z]);
						}
						else if ( z == m )
						{
							tmp[z] = f*cos(x[z]);
							f *= -sin(x[z]);
						}
						else
						{
							tmp[z] = f*sin(x[z]);
							f *= cos(x[z]);
						}
					}
					tmp[d-1] = f;

					dx[m] -= v_dot_pdt(grid_pts[k], tmp, d)/dist;
				}
			}
			for ( k = 0; k < d-1; k++ )
				dx[k] *= rate;

			// new value of grid_pts[j]
			f = 1;
			for ( k = 0; k < d-1; k++ )
			{
				grid_pts[j][k] = f * sin(x[k]+dx[k]);
				f *= cos(x[k]+dx[k]);
			}
			grid_pts[j][d-1] = f;
		}
	}
	delete[] x;
	delete[] dx;
	delete[] tmp;

	delete[] grid_pts;
	return gridpts_buf;
}


/* used as a subroutine for computing robust coresets; only slightly different from getCoreSet4 */
/* return: the coreset */
glPointSet * glPointSet::peelOneLevel( int grids, double eps, double * gridpts )
{
	int k, i;

	unsigned char * flag;
	flag = ( unsigned char * ) malloc( sizeof ( unsigned char ) * num );
	assert( flag != NULL );
	memset( flag, 0 , sizeof ( unsigned char ) * num );

	// compute coreset
	glPointSet *coreset = getCoreSet4( grids, eps, flag, gridpts );

	glReal * buffer = create_buffer( dim );

	// remove points of coreset from this point set
	for ( k = 1, i = 0; i < num; i++ )
		if ( flag[ i ] == 0 )
		{
			getPoint( i + 1, buffer );
			setPoint( k, buffer );
			k ++;
		}
	// change point set size
	num = k-1;

	/* wrap up */
	release_buffer( buffer );
	free( flag );

	return coreset;
}


/* compute robust coreset up to level k */
glPointSet * glPointSet::getRobustCoreset( int grids, double eps, int k )
{
	int i;
	glPointSet p, *coreset;

	p.init( dim, num, pnt_set );
	coreset = NULL;

	double *grid_pts = uniform_sph( dim, grids );

	// iterative peeling, for k times
	for ( i = 0; i < k; i++ )
	{
		if ( coreset == NULL )
			coreset = p.peelOneLevel( grids, eps, grid_pts );
		else
			coreset->merge( p.peelOneLevel( grids, eps, grid_pts ) );
	}

	p.dump();
	return coreset;
}

bool glPointSet::merge( glPointSet *p )
{
	if ( dim != p->dim || p == NULL )
		return false;

	glReal *old_pnt_set = pnt_set;
	int old_num = num;
	delete[] pnt_array;
	
	init( dim, num + p->num );
	memcpy( pnt_set, old_pnt_set, sizeof( glReal ) * dim * old_num );
	memcpy( pnt_set + dim * old_num, p->pnt_set, sizeof( glReal ) * dim * p->num );

	release_buffer( old_pnt_set);
	p->dump();

	return true;
}

glReal glPointSet::getWidth( glVector *v, int k, glReal *w )
{
        glReal *lmin, *lmax, l, * buffer;
	int i,j;

	lmin = create_buffer( k );
	lmax = create_buffer( k);
        buffer = create_buffer( dim );

	for ( i = 0; i < k; i++ )
	{
		lmin[i] = 1e10;
		lmax[i] = -(1e10);
	}

        for ( i = 1; i <= num; i++ )
        {
                getPoint( i, buffer );
                l = - v->dotProduct( buffer ) / v->getLength();

		j = k-1;
		while ( j >= 0 && lmax[j] < l )
		{
			if ( j < k-1 )	lmax[j+1] = lmax[j];
			j--;
		}
		if ( j < k-1 )	lmax[j+1] = l;

		j = k-1;
		while ( j >= 0 && lmin[j] > l )
		{
			if ( j < k-1 ) lmin[j+1] = lmin[j];
			j--;
		}
		if ( j < k-1 )	lmin[j+1] = l;
        }

	for ( i = 0; i < k; i++ )
		w[i] = lmax[i]-lmin[i];

	release_buffer( lmin );
	release_buffer( lmax );
        release_buffer( buffer );

	return w[0];
}
