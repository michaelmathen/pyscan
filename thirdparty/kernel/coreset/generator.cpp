/*********************************************************\
 *
 * generator.cpp 
 *     - generate various point distributions
 *
 * Hai Yu (fishhai@cs.duke.edu), Nov 2003
 *
\*********************************************************/

#include <stdlib.h>
#include <math.h>
#include <myinc.h>
#include <string.h>
#include <stdio.h>
#include <generator.h>

#define PRIME	100003
#define HALF	50001

double ** cylindrical_surface( int n, int d, double height )
{
	double r = 1;				// radius assumed to be 1
	double **pts;
	int i,j;
	double len;

	pts = newbuf( n, d );

	for ( i = 0; i < n; i ++ )
	{
		// generate first d-1 coordinates
		len = 0;
		while ( EQUAL( len, 0 ) )
		{
			for ( j = 0; j < d-1; j ++ )
				pts[i][j] = (rand()%PRIME - HALF)*1.0 / HALF;
			len = v_len( pts[i], d-1 );
		}
		for ( j = 0; j < d-1; j ++ )
			pts[i][j] /= len;

		// generate the last coordinate
		pts[i][d-1] = (rand()%PRIME - HALF)*0.5* height / HALF;
	}

	return pts;
}

double ** unit_sphere( int n, int d )
{
	double **pts;
	int i,j;
	double len;

	pts = newbuf( n, d );

	for ( i = 0; i < n; i ++ )
	{
		len = 0;
		while ( EQUAL( len, 0 ) )
		{
			for ( j = 0; j < d; j ++ )
				pts[i][j] = (rand()%PRIME - HALF)*1.0 / HALF;
			len = v_len( pts[i], d );
		}
		for ( j = 0; j < d; j ++ )
			pts[i][j] /= len;
	}

	return pts;
}

double ** concentric( int n1, int n2, int n3, int d )
{
	double ** pts;
	int i,j;
	double len;

	pts = newbuf( n1+n2+n3, d );

	for ( i = 0; i < n1; i ++ )
        {
                len = 0;
                while ( EQUAL( len, 0 ) )
                {
                        for ( j = 0; j < d; j ++ )
                                pts[i][j] = (rand()%PRIME - HALF)*1.0 / HALF;
                        len = v_len( pts[i], d );
                }
                for ( j = 0; j < d; j ++ )
                        pts[i][j] /= len;
        }

	for ( i = n1; i < n1+n2; i ++ )
        {
                len = 0;
                while ( EQUAL( len, 0 ) )
                {
                        for ( j = 0; j < d; j ++ )
                                pts[i][j] = (rand()%PRIME - HALF)*1.0 / HALF;
                        len = v_len( pts[i], d )/2.0;
                }
                for ( j = 0; j < d; j ++ )
                        pts[i][j] /= len;
        }

	for ( i = n1+n2; i < n1+n2+n3; i ++ )
        {
                len = 0;
                while ( EQUAL( len, 0 ) )
                {
                        for ( j = 0; j < d; j ++ )
                                pts[i][j] = (rand()%PRIME - HALF)*1.0 / HALF;
                        len = v_len( pts[i], d )/3.0;
                }
                for ( j = 0; j < d; j ++ )
                        pts[i][j] /= len;
        }

	return pts;
}

double ** annulus( int n, int d, double ar, double ir)
{
	double **pts;
	int i,j;
	double len, r;

	pts = newbuf( n, d );

	for ( i = 0; i < n; i ++ )
	{
		len = 0;
		while ( EQUAL( len, 0 ) )
		{
			for ( j = 0; j < d; j ++ )
				pts[i][j] = (rand()%PRIME - HALF)*1.0 / HALF;
			len = v_len( pts[i], d );
		}
		if ( i == 0 ) r = ir;
		else if ( i == 1 ) r = ar; 
		else r = (rand()%PRIME)*(ar-ir) / PRIME + ir;

		for ( j = 0; j < d; j ++ )
			pts[i][j] = pts[i][j]*r/len;
	}

	return pts;
}


double ** clustered( int n, int d )
{
	double **pts;
	double r;
	int i,j,k,l;

	double *center = new double[d];
	pts = newbuf( n, d );

	k = 0;
	for ( i = 1; i <= 20; i ++ )
	{
		for ( j = 0; j < d; j ++ )
			center[j] = (rand()%PRIME)*1.0 / PRIME;
		r  = (rand()%PRIME + 1.0) / PRIME / 5;

		l = 1;
		while ( l <= n/20 )
		{
			for ( j = 0; j < d; j++ )
				pts[k][j] = (rand()%PRIME  - HALF) * r / HALF;
			if ( v_len( pts[k], d ) > r ) continue;
			for ( j = 0; j < d; j++ )
				pts[k][j] += center[j];
			l++;
			k++;
		}
	}

	delete[] center;
	return pts;
}


double ** box3d_boundary( int n, double h, double l )   // width = 1
{
	double **pts;
	double b1,b2,b3;
	int face;

	pts = newbuf( n, 3 );

	for ( int i = 0; i < n; i ++ )
	{
		b1 = (rand() % PRIME)*(l*h+l+h) / PRIME;
		if ( b1 < l*h )
			b1 = l, b2 = h, b3 = 1, face = 1;
		else if ( b1 < l*h+h)
			b1 = h, b2 = 1, b3 = l, face = 2;
		else
			b1 = 1, b2 = l, b3 = h, face = 3;

		pts[i][ face%3 ] = (rand() % PRIME - HALF)*b1 / HALF;
		pts[i][ (face+1)%3 ] = (rand() % PRIME - HALF)*b2 / HALF;
		pts[i][ (face+2)%3 ] = b3 * ( (rand()%2)*2-1 );
	}

	return pts;
}

double ** ellipsoid3d( int n, double a, double b, double c )
{
	double **pts;
	double r;


	pts = newbuf( n, 3 );

	for ( int i = 0; i < n; i ++ )
	{
		pts[i][0] = (rand() % PRIME - HALF)*a / HALF + (rand() % PRIME -HALF)*0.1/HALF;
		r = sqrt( 1-SQR( pts[i][0]/a ) ) * b ;
		pts[i][1] = (rand() % PRIME - HALF)*r / HALF + (rand() % PRIME - HALF)*0.1/HALF;
		r = sqrt( 1-SQR( pts[i][0]/a )-SQR( pts[i][1]/b ) ) * c;
		pts[i][2] = ((rand() % 2)* 2 -1) * r + (rand() % PRIME - HALF) * 0.1/HALF;
	}

	return pts;
}

void dump( double** pts )
{
	delete[] pts[0];
	delete[] pts;
}

double** newbuf( int n, int d )
{
	double* pts_buf = new double[n*d];
	double** pts = new double*[n];
	for ( int i = 0; i < n; i ++ )
		pts[i] = &(pts_buf[i*d]);
	return pts;
}

double ** read_plyfile( char* filename, int &n )
{
	int i, j;
	float temp;
	double cd[3];
	char * buffer = new char[800];

	FILE * data = fopen( filename, "r" );
	if (data == NULL)
	{
		printf("error opening input file.\n");
		return NULL;
	}

	buffer[0] = 0;
	while ( strcmp( buffer, "element vertex" ) != 0 )
	{
		fgets(buffer, 800, data);
		buffer[14] = 0;
	}
	sscanf( buffer + 15, "%d", &n);
	double **pts = newbuf( n, 3 );

	while ( strcmp( buffer, "end_header") != 0 )
	{
		fgets(buffer, 800, data);
		buffer[10] = 0;
	}

	for ( i = 0; i < n; i++ )
	{
		fgets(buffer, 800, data);

		sscanf( buffer, "%f", &temp );
		pts[i][0] = temp;
		j = 0;
		while ( buffer[j] != 32 ) j++;
		j++;
		sscanf( buffer+j, "%f", &temp );
		pts[i][1] = temp;
		while ( buffer[j] != 32 ) j++;
		j++;
		sscanf( buffer+j, "%f", &temp );
		pts[i][2] = temp;
	}

	fclose( data );
	delete[] buffer;
	return pts;
}

