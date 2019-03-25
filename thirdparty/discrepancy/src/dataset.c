/*****************************
 * dataset.c
 * Functions on a dataset
 *
 * Jeff Phillips
 * jeffp@cs.duke.edu
 * 8.02.05
 *****************************/

#include <stdio.h>
#include "dataset.h"
#include "sort.h"
#include "mtint.h"

 int parent (int i) {return (i+1)/2 -1;}
 int left   (int i) {return 2*(i+1) -1;}
 int right  (int i) {return 2*(i+1);}

/**
 *  Reads in a file from a file.  Builds data.
 * 
 * @param ds    Structure to house data
 * @param ifp   Input file pointer
 * @return      0 if ok, 1 if file was not large enough
 */
int read_ds (data_set *ds, FILE *ifp) {
  int i,ok;
  char header[128];

  (void)fscanf(ifp, "%s", header);
  //fscanf(ifp, "%d", &length);
  //if (length < ds->npts) return 1;
  ds->data = (data_point*) malloc (sizeof(data_point) * ds->npts);

  for (i=0; i<ds->npts; ++i) {
    ok = fscanf(ifp, "%lf, %lf, %lf, %lf",
		&ds->data[i].x, &ds->data[i].y, &ds->data[i].blue, &ds->data[i].red);
    if (!ok) return 1;
    ds->rpts += ds->data[i].red;
    ds->bpts += ds->data[i].blue;
  }
  return 0;
}

/** 
 *  Writes a file from dataset.
 * 
 * @param ds      Dataset to output to a file.
 * @param ofp     Output file pointer
 */
void write_ds (data_set *ds, FILE *ofp) {
  int i; 
  (void)fprintf(ofp, "x,y,pop,count\n\n");
  //fprintf(ofp, "%d\n\n", ds->npts);
  
  for (i=0; i<ds->npts; ++i) 
    fprintf (ofp, "%f, %f, %f, %f\n", ds->data[i].x, ds->data[i].y, ds->data[i].blue, ds->data[i].red);
}


/**
 * Creates and loads a grid with data points.
 * 
 * @param ds          Dataset to use to fill grid.
 * @param grid_size   Size of each dimension of grid.
 * @param red         Pointer to unallocated red grid.
 * @param blue        Pointer to unallocated blue grid.
 */
void load_grid_ds (data_set *ds, int grid_size, double ***red, double ***blue,
		   double *minx, double *dx, double *miny, double *dy) {
  int x,y,k;

  double maxx=ds->data[0].x;
  double maxy=ds->data[0].y;

  *minx = ds->data[0].x;
  *miny = ds->data[0].y;

  *red = (double**) malloc (sizeof(double*) * grid_size);
  *blue = (double**) malloc (sizeof(double*) * grid_size);
  for (x=0; x<grid_size; ++x) {
    (*red)[x] = (double*) malloc (sizeof(double) * grid_size);
    (*blue)[x] = (double*) malloc (sizeof(double) * grid_size);
    for (y=0; y<grid_size; ++y) {
      (*red)[x][y] = 0.0;
      (*blue)[x][y] = 0.0;
    }
  }

  for (k=0; k<ds->npts; ++k) {
    if (ds->data[k].x < *minx) *minx = ds->data[k].x;
    if (ds->data[k].x > maxx) maxx = ds->data[k].x;
    if (ds->data[k].y < *miny) *miny = ds->data[k].y;
    if (ds->data[k].y > maxy) maxy = ds->data[k].y;
  }

  *dx = (double)(maxx*1.00001-(*minx))/grid_size;
  *dy = (double)(maxy*1.00001-(*miny))/grid_size;

  for (k=0; k<ds->npts; ++k) {
    x = (int)((ds->data[k].x - (*minx))/(*dx));
    y = (int)((ds->data[k].y - (*miny))/(*dy));
    (*red)[x][y] += ds->data[k].red;
    (*blue)[x][y] += ds->data[k].blue;
  }
}



/**
 *  Builds a random dataset from 2 Poisson distributions
 * 
 * @param ds       Structure to house data
 * @param n_pts    Number of points in set
 * @param mean1    Mean of first distribution
 * @param stdev1   Standeard Deviation of first distribution
 * @param mean2    Mean of second distribution
 * @param stdev2   Standard Deviation of second distribution
 * @param x_min    [0,1], < x_max - min x coord of second set
 * @param x_max    [0,1], > x_min - max x coord of second set
 * @param y_min    [0,1], < y_max - min y coord of second set
 * @param y_max    [0,1], > y_min - max y coord of second set
 */
void random_Gaussian2_ds (data_set *ds, int n_pts,
			  double mean1, double stdev1, double mean2, double stdev2, 
			  double x_min, double x_max, double y_min, double y_max) {
  int i;
  double mean = 100.0;
  double stdev = 50.0;
  ds->npts = n_pts;

  sgenrand_auto();
  ds->data = (data_point*) malloc (sizeof(data_point) * ds->npts);
  for (i=0; i<n_pts; ++i) {
    // position
    ds->data[i].x = genuniform();
    ds->data[i].y = genuniform();
    // blue
    ds->data[i].blue = mean + stdev*gennormal();
    if (ds->data[i].blue < 0) ds->data[i].blue = 0;
    ds->bpts += ds->data[i].blue;
    // red
    if (ds->data[i].x > x_min && ds->data[i].x < x_max && 
	ds->data[i].y > y_min && ds->data[i].y < y_max)
      ds->data[i].red = mean2 + stdev2*gennormal();
    else
      ds->data[i].red = mean1 + stdev1*gennormal();
    ds->rpts += (int)ds->data[i].red;
  }

}


/**
 *  Randomly builds a dataset -- Gaussian style
 * 
 * @param ds     Structure to house data
 * @param n_pts  Number of data points
 */
void random_Gaussian_ds (data_set *ds, int n_pts) {
  random_Gaussian2_ds(ds, n_pts, 40, 10, 40, 10, 0, 0, 0, 0);
}


void random_Poisson_file_ds(data_set *ds, int n_pts, double frac1, double frac2, 
			    double x_min, double x_max, double y_min, double y_max, 
			    FILE *ifp) {
  int SIZE = 61291;
  int i,j,r,b;
  double x,y;
  char header[128];
  double X[SIZE];
  double Y[SIZE];
  ds->npts = n_pts;
  sgenrand_auto();
  ds->data = (data_point*) malloc (sizeof(data_point) * ds->npts);
  
  // read in (x,y) from file.
    (void)fscanf(ifp, "%s", header);

  for (i=0; i<SIZE; ++i) {
    (void)fscanf(ifp, "%lf, %lf, %d, %d", &x, &y, &r, &b);
    X[i] = (x-genuniform())/256.0;
    Y[i] = (y-genuniform())/256.0;
  }
  
  for (i=0; i<n_pts; ++i) {
    j = (int)(genuniform()*SIZE);   
    ds->data[i].x = X[j];
    ds->data[i].y = Y[j];
    ds->data[i].blue = (int)exp(genuniform()*6);
    if (ds->data[i].blue < 0) ds->data[i].blue = 0;
    ds->bpts += ds->data[i].blue;
    ds->data[i].red = 0;
    if (ds->data[i].x > x_min && ds->data[i].x < x_max && 
	ds->data[i].y > y_min && ds->data[i].y < y_max) 
      for (j=0; j<ds->data[i].blue; ++j) {
	if (genuniform() < frac2)
	  ds->data[i].red += 1;
      }
    else
      for (j=0; j<ds->data[i].blue; ++j)
	if (genuniform() < frac1)
	  ds->data[i].red += 1;
    ds->rpts += (int)ds->data[i].red;
  }   
}


/**
 *  Builds a random dataset from 2 Poisson distributions
 * 
 * @param ds       Structure to house data
 * @param n_pts    Number of points in set
 * @param frac1    Fraction of first set to be red
 * @param frac2    Fraction of second set to be red
 * @param x_min    [0,1], < x_max - min x coord of second set
 * @param x_max    [0,1], > x_min - max x coord of second set
 * @param y_min    [0,1], < y_max - min y coord of second set
 * @param y_max    [0,1], > y_min - max y coord of second set
 */
void random_Poisson2_ds (data_set *ds, int n_pts, double frac1, double frac2,
			 double x_min, double x_max, double y_min, double y_max) {
  int i,j;
  //  double mean = 100.0;
  //  double stdev = 50.0;
  ds->npts = n_pts;

  sgenrand_auto();
  ds->data = (data_point*) malloc (sizeof(data_point) * ds->npts);
  for (i=0; i<n_pts; ++i) {
    ds->data[i].x = genuniform();
    ds->data[i].y = genuniform();
    ds->data[i].blue = (int)exp(genuniform()*6);
    if (ds->data[i].blue < 0) ds->data[i].blue = 0;
    ds->bpts += ds->data[i].blue;
    ds->data[i].red = 0;
    if (ds->data[i].x > x_min && ds->data[i].x < x_max && 
	ds->data[i].y > y_min && ds->data[i].y < y_max) 
      for (j=0; j<ds->data[i].blue; ++j) {
	if (genuniform() < frac2)
	  ds->data[i].red += 1;
      }
    else
      for (j=0; j<ds->data[i].blue; ++j)
	if (genuniform() < frac1)
	  ds->data[i].red += 1;
    ds->rpts += (int)ds->data[i].red;
  }
}


/**
 *  Randomly builds a dataset - Poisson style
 * 
 * @param ds     Structure to house data
 * @param n_pts  Number of data points
 */
void random_Poisson_ds (data_set *ds, int n_pts) {
  random_Poisson2_ds(ds, n_pts, .05, .05, 0, 0, 0, 0);
}



/**
 *  Randomly builds a dataset - Bernoulli style
 * 
 * @param ds     Structure to house data
 * @param n_pts  Number of data points
 */
void random_Bernoulli_ds (data_set *ds, int n_pts) {
  int i,j;
  double frac = .4;
  int b_mean = 15;
  double stdev = 3;
  ds->npts = n_pts;

  sgenrand_auto();
  ds->data = (data_point*) malloc (sizeof(data_point) * ds->npts);
  for (i=0; i<ds->npts; ++i) {
    ds->data[i].x = genuniform();
    ds->data[i].y = genuniform();
    ds->data[i].blue = (int)(b_mean + stdev*gennormal());
    if (ds->data[i].blue < 0) ds->data[i].blue = 0;
    ds->bpts += ds->data[i].blue;
    ds->data[i].red = 0;
    for (j=0; j<ds->data[i].blue; ++j)
      if (genuniform() < frac)
	ds->data[i].red += 1;
    ds->rpts += (int)ds->data[i].red;
  }
}

/**
 *  Randomly builds a dataset - gamma style
 *  actually Poisson - cause I am lazy...
 * 
 * @param ds     Structure to house data
 * @param n_pts  Number of data points
 */
void random_gamma_ds (data_set *ds, int n_pts) {
  int i,j;
  double frac = .05;
  int b_mean = 15;
  double stdev = 3;
  ds->npts = n_pts;

  sgenrand_auto();
  ds->data = (data_point*) malloc (sizeof(data_point) * ds->npts);
  for (i=0; i<ds->npts; ++i) {
    ds->data[i].x = genuniform();
    ds->data[i].y = genuniform();
    ds->data[i].blue = (b_mean + stdev*gennormal());
    if (ds->data[i].blue < 0) ds->data[i].blue = 0;
    ds->bpts += ds->data[i].blue;
    ds->data[i].red = 0;
    for (j=0; j<ds->data[i].blue; ++j)
      ds->data[i].red += genuniform()*frac;
    ds->rpts += ds->data[i].red;
  }
}


/**
 * Sets the red and blue weight of the points to 
 *       represent a specific hyperplane.
 *
 * @param ds    Data_set
 * @param rw    red weight
 * @param bw    blue weight
 */
void set_red_blue_weight_ds (data_set *ds, double rw, double bw) {
  ds->red_weight = rw;
  ds->blue_weight = bw;
}

/**
 * Creates the heap structure, does not fill.
 * 
 * @param ds    data set
 */
void create_heap_ds (data_set *ds) {
  int i;
  int base_size=1;
  int size;
  double *xpts = (double*) malloc (sizeof(double)*ds->npts);
  double *ypts = (double*) malloc (sizeof(double)*ds->npts);
  while (base_size < ds->npts) base_size*=2;
  size = base_size+ds->npts+1;
  ds->iheap = (ival*) malloc (sizeof(ival) * size);
  ds->xorder = (int*) malloc (sizeof(int) * ds->npts);
  ds->yorder = (int*) malloc (sizeof(int) * ds->npts);

  // sort points
  for (i=0; i<ds->npts; ++i) {
    ds->xorder[i] = i;  
    ds->yorder[i] = i;
    xpts[i] = ds->data[i].x;
    ypts[i] = ds->data[i].y;
  }
  idquicksort (xpts, ds->xorder, 0, ds->npts-1);
  idquicksort (ypts, ds->yorder, 0, ds->npts-1);
  free (xpts);
  free (ypts);
  for (i=0; i<ds->npts; ++i) {
    ds->data[ds->xorder[i]].xorder = i;
    ds->data[ds->yorder[i]].yorder = i;
  }

  // set interval boundaries
  ds->iheap[size-1].rB = ds->data[ds->yorder[ds->npts-1]].y;
  ds->iheap[base_size].lB = ds->data[ds->yorder[0]].y;
  for (i=0; i<ds->npts-1; ++i) {
    double mid = (ds->data[ds->yorder[i]].y + ds->data[ds->yorder[i+1]].y)/2.0;
    ds->iheap[base_size+i].rB = mid;
    ds->iheap[base_size+i+1].lB = mid;
  }
  for (i=parent(size-1); i>=0; --i) {
    ds->iheap[i].lB = ds->iheap[left(i)].lB;
    ds->iheap[i].rB = ds->iheap[right(i)].rB;
  }

  // zero out the rest of it
  for (i=0; i<size; ++i) {
    clear_ival(&(ds->iheap[i]));
  }			     
}


/**
 * Creates the heap structure for on a grid, does not fill.
 * 
 * @param ds          data set
 * @param gsize       size of data in heap (number of leaves)
 */
void create_grid_heap_ds (data_set *ds, int gsize) {
  int i;
  int base_size=1;
  int size;
  while (base_size < gsize) base_size*=2;
  size = base_size+gsize-1;
  ds->iheap = (ival*) malloc (sizeof(ival) * size);

  // set interval boundaries
  ds->iheap[base_size+gsize-1].rB = gsize;
  ds->iheap[base_size].lB = 0;
  for (i=0; i<base_size-1; ++i) {
    ds->iheap[base_size+i].rB = i+1;
    ds->iheap[base_size+i+1].lB = i+1;
  }
  for (i=parent(size-1); i>=0; --i) {
    ds->iheap[i].lB = ds->iheap[left(i)].lB;
    ds->iheap[i].rB = ds->iheap[right(i)].rB;
  }

  // zero out the rest of it
  for (i=0; i<size; ++i) {
    clear_ival(&(ds->iheap[i]));
  }			     
}


/**
 * Computes the maximum discrepancy rectangle
 * 
 * @param ds   data set
 * @param xmin  reference to min x of range
 * @param xmax  reference to max x of range
 * @param ymin  reference to min y of range
 * @param ymax  reference to max y of range
 * @param nr    amount of red points in range
 * @param nb    amount of blue points in range
 * @return      Maximum discrepancy
 */
double max_discrepancy_ds (data_set *ds, 
			   double *xmin, double *xmax, double *ymin, double *ymax,
			   double *nr, double *nb) {
  int x1;
  int x2;
  double maxd = 0;
  int i;
  int base_size=1;
  int size;
  *xmin=0; *xmax=0; *ymin=0; *ymax=0;

  while (base_size < ds->npts) base_size*=2;
  size = base_size+ds->npts+1;

  for (x1=0; x1<ds->npts; ++x1) {

    // reset heap
    for (x2=0; x2<size; ++x2)
      clear_ival(&(ds->iheap[x2]));

    // rebuild heap
    for (x2=x1; x2<ds->npts; ++x2) {

      i = ds->data[ds->xorder[x2]].yorder;

      add_point_ival(&(ds->iheap[base_size+i]), 
		     ds->data[ds->xorder[x2]].y,
		     ds->data[ds->xorder[x2]].red,
		     ds->data[ds->xorder[x2]].blue,
		     ds->red_weight, ds->blue_weight);

      for (i = parent(base_size+i); i>=0; i = parent(i)) {
	merge_ival (&(ds->iheap[i]), 
		    &(ds->iheap[left(i)]), 
		    &(ds->iheap[right(i)]), ds->rpts, ds->bpts);
	//printf ("(%d %2.2f) ", i, ds->iheap[i].m.d);
      }

      //printf ("|%d %d - max discrepancy %f\n", x1, x2, ds->iheap[0].m.d);
      if (ds->iheap[0].m.d > maxd) {
      	maxd = ds->iheap[0].m.d;
      	*xmin = ds->data[ds->xorder[x1]].x;
      	*xmax = ds->data[ds->xorder[x2]].x;
      	*ymin = ds->iheap[0].m.l;
      	*ymax = ds->iheap[0].m.r;
      	*nr = ds->iheap[0].nr;
      	*nb = ds->iheap[0].nb;
      }
    }
  }

  return maxd;
}


/**
 * Computes the maximum discrepancy rectangle on a grid = [gsize]x[gsize]
 * 
 * @param ds     data set
 * @param gsize  size of data
 * @param red    red points grid
 * @param blue   blue points grid
 * @param xmin   reference to min x of range
 * @param xmax   reference to max x of range
 * @param ymin   reference to min y of range
 * @param ymax   reference to max y of range
 * @return       Maximum discrepancy
 */
double max_grid_discrepancy_ds (data_set *ds, int gsize, double **red, double **blue,
				double *xmin, double *xmax, double *ymin, double *ymax) {
  int x1;
  int x2;
  int y;
  double maxd = 0;
  int i;
  int base_size=1;
  int size;
  *xmin=0; *xmax=0; *ymin=0; *ymax=0;

  while (base_size < gsize) base_size*=2;
  size = base_size+gsize-1;

  for (x1=0; x1<gsize; ++x1) {

    // reset heap
    for (x2=0; x2<size; ++x2)
      clear_ival(&(ds->iheap[x2]));

    // rebuild heap
    for (x2=x1; x2<gsize; ++x2) {

      //i = ds->data[ds->xorder[x2]].yorder;

      for (y=0; y<gsize; ++y) 
	add_point_to_ival(&(ds->iheap[base_size+y]), y, red[x2][y], blue[x2][y],
			  ds->red_weight, ds->blue_weight);

      for (i = parent(size-1); i>=0; --i)
	merge_ival (&(ds->iheap[i]), 
		    &(ds->iheap[left(i)]), 
		    &(ds->iheap[right(i)]), ds->rpts, ds->bpts);

      if (ds->iheap[0].m.d > maxd) {
	maxd = ds->iheap[0].m.d;
	*xmin = x1;
	*xmax = x2;
	*ymin = ds->iheap[0].m.l;
	*ymax = ds->iheap[0].m.r;
      }
    }
  }

  return maxd;
}
