/*********************************************************\
 *
 * generator.h 
 *     - header file for generator.cpp
 *       generate various point distributions
 *
 * Hai Yu (fishhai@cs.duke.edu), Nov 2003
 *
\*********************************************************/

#ifndef __GENERATOR__H
#define __GENERATOR__H

/* generate points from a cylindrical surface of radius 1 */
double ** cylindrical_surface( int n, 		/* number of points */
                               int d, 		/* dimension */
                               double h ); 	/* height */

/* generate points from a unit sphere of radius 1 */
double ** unit_sphere( int n, 			/* number of points */
                       int d );			/* dimension */

/* generate clustered points */
double ** clustered( int n, 			/* number of points */
                     int d );			/* dimension */

/* generate points from an annulus */
double ** annulus( int n,			/* number of points */ 
                   int d, 			/* dimension */
                   double ar, 			/* outer radius */
                   double ir );			/* inner radius */

/* generate points from the surface of a 3d box of width 1 */
double ** box3d_boundary( int n, 		/* number of points */
                          double h, 		/* height of box */
                          double l );		/* length of box */

/* generate points from a 3d ellipsoid surface */
double ** ellipsoid3d( int n, 			/* number of points */
                       double a, double b, double c ); /* radii */

/* generate points from three concentric spheres */
double ** concentric( int n1, 			/* # of points on sphere of radius 1 */
                      int n2, 			/* # of points on sphere of radius 1/2 */
                      int n3,			/* # of points on sphere of radius 1/3 */ 
                      int d );			/* dimension */

/* read .ply format file */
double ** read_plyfile( char* filename, 	/* file name */
                        int &n );		/* # of points returned */

/* release array of points */
void   dump( double ** pts );			/* release pts[0] and pts */

/* create array of points */
double ** newbuf( int n, 			/* number of points */
                  int d );			/* dimension */

#else  /* __GENERATOR__H */
#error header file generator.h included twice
#endif /* __GENERATOR__H */

