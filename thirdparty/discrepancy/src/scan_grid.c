/************************
 * scan-grid.c
 * Finds the rectangle that maximizes some statistical discrepancy
 * 
 * Jeff M Phillips
 * jeffp@cs.duke.edu
 * 12.18.05
 ************************/

#include "scan_grid.h"

/**
 * Finds the rectange with the largest discrepancy by checking all rectangles
 * 
 * @param ds    Data set
 * @param f     Discrepancy function
 * @param xmin  reference to min x of range
 * @param xmax  reference to max x of range
 * @param ymin  reference to min y of range
 * @param ymax  reference to max y of range
 * @return      Maximum discrepancy
 */
double find_rect_grid (data_set *ds, d_fxn *f, int grid_size, double **red, double **blue,
		       double *xmin, double *xmax, double *ymin, double *ymax) {
  int x1, x2, y1, y2, y;
  double r,b;

  double disc, maxd = 0;
  *xmin=0; *xmax=0; *ymin=0; *ymax=0;

  double *r1 = (double*) malloc (sizeof (double) * grid_size);
  double *b1 = (double*) malloc (sizeof (double) * grid_size);

  for (x1=0; x1 < grid_size; ++x1) {   //(1)
    for (y=0; y<grid_size; ++y) {  
      r1[y] = 0.0;    
      b1[y] = 0.0;
    }
    for (x2=x1; x2 < grid_size; ++x2) {   //(2)
      if (x2==x1) {
	r1[0] = red[x2][0]; 
	b1[0] = blue[x2][0];
	for (y=1; y<grid_size; ++y) {
	  r1[y] = r1[y-1]+ red[x2][y];
	  b1[y] = b1[y-1]+blue[x2][y];
	}
      } 
      else {
	r=0; b=0;
	for (y=0; y<grid_size; ++y) {
	  r += red[x2][y];       r1[y] += r;
	  b += blue[x2][y];      b1[y] += b;
	}
      }
      
      for (y1=0; y1 < grid_size; ++y1)    //(3)
	for (y2=y1; y2 < grid_size; ++y2) {	//(4)  
	  if (y1!=0) {
	    r = r1[y2] - r1[y1-1];
	    b = b1[y2] - b1[y1-1];
	  }
	  else {
	    r = r1[y2];
	    b = b1[y2];
	  }

	  if (r>=BASE && b>=BASE && ds->rpts-r>=BASE && ds->bpts-b>=BASE &&
	      (double)r/ds->rpts > (double)b/ds->bpts) {
	    disc = f->f((double)r/ds->rpts, (double)b/ds->bpts);
	    if (disc > maxd)  {
	      maxd = disc;
	      *xmin = x1;
	      *xmax = x2;
	      *ymin = y1;
	      *ymax = y2;
	    }
	  }
	}
    }
  }
  free(r1);  
  free(b1);
  return maxd;
}
