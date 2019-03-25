/************************
 * stat-rect.c
 * Finds the rectangle that maximizes some statistical discrepancy
 * 
 * Jeff M Phillips
 * jeffp@cs.duke.edu
 * 8.08.05
 ************************/

#include "stat_rect.h"



//inline void print_usage() {
//  printf (">> stat-rect 'dist-type={P,G,B,g}' #error #num-pts <infile>\n");
//  exit(2);
//}
//
//int main (int argc, char* argv[]) {
//  int i;
//
//  double x,y;
//  double xmin, xmax, ymin, ymax, maxnr, maxnb;
//  double x_min, x_max, y_min, y_max;
//  FILE *ifp;
//
//  d_fxn D;
//  ad_fxn A;
//
//  char dist_type = 'P';
//  int n_pts = 100;
//  double eps = .01;
//
//  data_set ds;
//  double m, maxd=0;
//  clock_t c1, c2;
//
//  if (argc > 5)   print_usage();
//  if (argc >= 4)  n_pts = atoi(argv[3]); ds.npts = n_pts;
//  if (argc >= 3)  eps = atof(argv[2]);
//  if (argc >= 2)  dist_type = argv[1][0];
//
//  // read file
//  if (argc == 5) {
//    ifp = fopen(argv[4], "r");
//    if (!ifp || read_ds(&ds, ifp)) {
//      printf ("Bad infile: %s\n", argv[4]);
//      print_usage();
//    }
//  }
//
//  // generate random dataset + set function type
//  if (dist_type == 'P') {
//    if (argc < 5) random_Poisson_ds(&ds, ds.npts);
//    set_Poisson_fxn(&D);
//  }
//  else if (dist_type == 'G') {
//    if (argc < 5) random_Gaussian_ds(&ds, ds.npts);
//    set_Gaussian_fxn(&D);
//  }
//  else if (dist_type == 'B') {
//    if (argc < 5) random_Bernoulli_ds(&ds, ds.npts);
//    set_Bernoulli_fxn(&D);
//  }
//  else if (dist_type == 'g') {
//    if (argc < 5) random_gamma_ds(&ds, ds.npts);
//    set_gamma_fxn(&D);
//  }
//  else print_usage();
//
//
//  // find approximate function
//  c1 = clock();
//  //get_gridding_fxn(&D, &A, eps, ds.rpts, ds.bpts);
//  //dumb_gridding_fxn(&D, &A, eps, ds.rpts, ds.bpts);
//  smart_gridding_fxn(&D, &A, eps, ds.npts);
//  c2 = clock();
//  printf ("got gridding: %d  -- time: %f\n", A.n, (double)(c2-c1)/CLOCKS_PER_SEC);
//  printf ("red: %f   blue: %f\n", ds.rpts, ds.bpts);
//
//  // test some points possible rectangles
///*   for (x=.1; x<.9; x+=.1) */
///*     for (y=.1; y<x; y+=.1) */
///*       printf("(%2.1f %2.1f)  D:%2.3f  >=  A:%2.3f \n", */
///* 	     x, y, D.f(x,y), eval_approx_fxn(&A, x, y)); */
//
//  // compute max discrepancy
//  c1 = clock();
//  create_heap_ds(&ds);
//  for (i=0; i<A.n; ++i) {
//    set_red_blue_weight_ds(&ds, A.x[i]/ds.rpts, A.y[i]/ds.bpts);
//    m = max_discrepancy_ds(&ds, &x_min, &x_max, &y_min, &y_max, &maxnr, &maxnb) + A.z[i];
//    if (m > maxd) {
//      maxd = m;
//      xmin = x_min;
//      xmax = x_max;
//      ymin = y_min;
//      ymax = y_max;
//    }
//  }
//  c2 = clock();
//  printf ("maximum discrepancy:   %f -- time:%f   [%1.2f,%1.2f]x[%1.2f,%1.2f]\n",
//	  maxd,  (double)(c2-c1)/CLOCKS_PER_SEC, xmin, xmax, ymin, ymax);
//
//  // naively check max discrepancy
//  c1 = clock();
//  m = naive_find_rect(&ds, &D, &xmin, &xmax, &ymin, &ymax);
//  c2 = clock();
//  printf ("naive max disc:        %f -- time:%f   [%1.2f,%1.2f]x[%1.2f,%1.2f]\n",
//	  m,  (double)(c2-c1)/CLOCKS_PER_SEC, xmin, xmax, ymin, ymax);
//
//  // naive check max with approximate discrepancy function
///*   c1 = clock(); */
///*   m = naive_approx_find_rect(&ds, &A, &xmin, &xmax, &ymin, &ymax); */
///*   c2 = clock(); */
///*   printf ("naive approx max disc: %f -- time:%f   [%1.2f,%1.2f]x[%1.2f,%1.2f]\n",  */
///* 	  m,  (double)(c2-c1)/CLOCKS_PER_SEC, xmin, xmax, ymin, ymax); */
//
//
//  return 0;
//}
//



/**
 * Finds the rectange with the largest discrepancy by checking all rectangles
 *         Uses approximate discrepancy function 
 *
 * @param ds    Data set
 * @param f     Approximate discrepancy function
 * @param xmin  reference to min x of range
 * @param xmax  reference to max x of range
 * @param ymin  reference to min y of range
 * @param ymax  reference to max y of range
 * @return      Maximum discrepancy
 */
double naive_approx_find_rect (data_set *ds, ad_fxn *A, 
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
	  //	  printf ("|%d %d %d %2.4f|", x1, x2, y1, maxd);
	  for (y2=y1; y2 < ds->npts; ++y2) 
	    if (ds->data[ds->yorder[y2]].xorder >= ds->data[ds->xorder[x1]].xorder && 
		ds->data[ds->yorder[y2]].xorder <= ds->data[ds->xorder[x2]].xorder) {
	      r += ds->data[ds->yorder[y2]].red;
	      b += ds->data[ds->yorder[y2]].blue;
	      if (r>=BASE && b>=BASE && ds->rpts-BASE>=r && ds->bpts-BASE>=b &&
		  (double)r/ds->rpts > (double)b/ds->bpts) {
		disc = eval_approx_fxn(A, (double)r/ds->rpts, (double)b/ds->bpts);
		//		printf("[%2.4f - %d] ", disc, disc>maxd);
		if (disc > maxd)  {
		  maxd = disc;
		  *xmin = ds->data[ds->xorder[x1]].x;
		  *xmax = ds->data[ds->xorder[x2]].x;
		  *ymin = ds->data[ds->yorder[y1]].y;
		  *ymax = ds->data[ds->yorder[y2]].y;
		}
	      }
	    }
	  //	  printf ("\n");
	}
  
  return maxd;

}
