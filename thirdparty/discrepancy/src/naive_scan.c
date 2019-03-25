/************************
 * naive-scan.c
 * Finds the rectangle that maximizes some statistical discrepancy
 * 
 * Jeff M Phillips
 * jeffp@cs.duke.edu
 * 12.17.05
 ************************/

#include "dataset.h"
#include "function.h"

#include "naive_scan.h"


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
double naive_find_rect (data_set *ds, d_fxn *f, 
			double *xmin, double *xmax, double *ymin, double *ymax) {
  int x1, x2, y1, y2;
  double r,b;

  double disc, maxd = 0;
  *xmin=0; *xmax=0; *ymin=0; *ymax=0;

  for (x1=0; x1 < ds->npts; ++x1) 
    for (x2=x1+1; x2 < ds->npts; ++x2) 
      for (y1=0; y1 < ds->npts; ++y1) 
	if (ds->data[ds->yorder[y1]].xorder >= ds->data[ds->xorder[x1]].xorder && 
	    ds->data[ds->yorder[y1]].xorder <= ds->data[ds->xorder[x2]].xorder) {
	  r=0; b=0;
	  //printf ("|%d %d %d %2.4f|", x1, x2, y1, maxd);
	  for (y2=y1; y2 < ds->npts; ++y2) 
	    if (ds->data[ds->yorder[y2]].xorder >= ds->data[ds->xorder[x1]].xorder && 
		ds->data[ds->yorder[y2]].xorder <= ds->data[ds->xorder[x2]].xorder) {
	      r += ds->data[ds->yorder[y2]].red;
	      b += ds->data[ds->yorder[y2]].blue;
	      if (r>=BASE && b>=BASE && ds->rpts-BASE>=r && ds->bpts-BASE>=b &&
		  (double)r/ds->rpts > (double)b/ds->bpts) {
		disc = f->f((double)r/ds->rpts, (double)b/ds->bpts);
		//printf("[%2.4f - %d] ", disc, disc>maxd);
		if (disc > maxd)  {
		  maxd = disc;
		  *xmin = ds->data[ds->xorder[x1]].x;
		  *xmax = ds->data[ds->xorder[x2]].x;
		  *ymin = ds->data[ds->yorder[y1]].y;
		  *ymax = ds->data[ds->yorder[y2]].y;
		}
	      }
	    }
	  //printf ("\n");
	}
  
  return maxd;

}
