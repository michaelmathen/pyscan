/************************
 * scan-stat.c
 * Finds the rectangle that maximizes some statistical discrepancy
 * 
 * Jeff M Phillips
 * jeffp@cs.duke.edu
 * 12.17.05
 ************************/

#include "dataset.h"
#include "time.h"
#include "function.h"


//inline void print_usage() {
//  printf (">> scan-stat <infile> #eps #num-pts <outfile>\n");
//  exit(2);
//}
//
//int main (int argc, char* argv[]) {
//  int i;
//
//  double x,y;
//  double xmin, xmax, ymin, ymax;
//  double x_min, x_max, y_min, y_max;
//  FILE *ifp;
//  FILE *ofp;
//
//  d_fxn D;
//  ad_fxn A;
//
//  int n_pts = 100;
//  double eps = .01;
//
//  data_set ds;
//  double m, maxd=0, maxnr, maxnb;
//  clock_t c1, c2;
//
//  if (argc > 5 || argc < 3)   print_usage();
//  if (argc >= 5)  ofp = fopen(argv[4], "w");
//  else            ofp = fopen("regions.csv", "w");
//  if (!ofp) {printf("bad output file: %s\n", argv[4]); print_usage();}
//  if (argc >= 4)  n_pts = atoi(argv[3]);    ds.npts = n_pts;
//  if (argc >= 3)  eps = atof(argv[2]);
//  if (argc >= 2)  ifp = fopen(argv[1], "r");
//
//  // read file
//  if (!ifp || read_ds(&ds, ifp)) {
//    printf ("bad input file: %s\n", argv[1]);
//    print_usage();
//  }
//  fclose(ifp);
//
//  set_Poisson_fxn(&D);
//
//  // find approximate function
//  c1 = clock();
//  smart_gridding_fxn(&D, &A, eps, ds.npts);
//  c2 = clock();
//  printf ("got gridding: %d  -- time: %f\n", A.n, (double)(c2-c1)/CLOCKS_PER_SEC);
//  printf ("red: %f   blue: %f\n", ds.rpts, ds.bpts);
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
//
//  // write to file
//  fprintf(ofp, "score,x_lo,y_lo,x_hi,y_hi,count,pop,total_count,total_pop,p_value\n\n");
//  fprintf(ofp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
//	  m, xmin, ymin, xmax, ymax, 1.0, 2.0, ds.rpts, ds.bpts, 1.0,
//	  (double)(c2-c1)/CLOCKS_PER_SEC);
//  fclose(ofp);
//
//
//  return 0;
//}
//
//
