/************************
 * dataset.h
 * Storage for data points
 * 
 * Jeff Phillips
 * jeffp@cs.duke.edu
 * 8.02.05
 ***********************/

#ifndef DATASET__H
#define DATASET__H

//#define RED 1
//#define BLUE 0

#include "interval.h"
#include "header.h"

struct data_point_t {
  double x;
  double y;
  double red;
  double blue;
  int xorder;
  int yorder;
};
typedef struct data_point_t data_point;

struct data_set_t {
  data_point* data;
  int npts;//, nrpts, nbpts;
  double rpts, bpts;

  int* xorder;
  int* yorder;

  ival* iheap;

  double red_weight;
  double blue_weight;
};
typedef struct data_set_t data_set;

///////// FUNCTIONS ///////////////////
int read_ds (data_set *ds, FILE *ifp);
void write_ds (data_set *ds, FILE *ofp);
void load_grid_ds (data_set *ds, int grid_size, double ***red, double ***blue,
		   double *minx, double *dx, double *miny, double *dy);
void random_Poisson_ds (data_set *ds, int n_pts);
void random_Poisson2_ds (data_set *ds, int n_pts, double frac1, double frac2,
			 double x_min, double x_max, double y_min, double y_max);
void random_Gaussian_ds (data_set *ds, int n_pts);
void random_Gaussian2_ds (data_set *ds, int n_pts, double mean1, double stdev1, 
			  double mean2, double stdev2, 
			  double x_min, double x_max, double y_min, double y_max);
void random_Bernoulli_ds (data_set *ds, int n_pts);
void random_gamma_ds (data_set *ds, int n_pts);
void set_red_blue_weight_ds (data_set *ds, double rw, double bw);
void create_heap_ds (data_set *ds);
void create_grid_heap_ds (data_set *ds, int grid_size);
double max_discrepancy_ds (data_set *ds, 
			   double *xmin, double *xmax, double *ymin, double *ymax,
			   double *nr, double *nb); 
double max_grid_discrepancy_ds (data_set *ds, int grid_size, double **red, double **blue,
				double *xmin, double *xmax, double *ymin, double *ymax); 


#endif
