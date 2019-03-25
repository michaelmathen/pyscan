/************************
 * find-rect.c
 * Finds the rectangle that maximizes discrepancy
 * 
 * Jeff M Phillips
 * jeffp@cs.duke.edu
 * 8.03.05
 ************************/

#include "dataset.h"
#include "time.h"


//inline void print_usage() {
//  printf (">> find-rect #num_pts\n");
//  exit(2);
//}
//
//int main (int argc, char* argv[]) {
//
//  FILE *ifp;
//
//  double xmin, xmax, ymin, ymax;
//  double red_weight = 1.0;
//  double blue_weight = -1.0;
//
//  int n_pts = 100;
//
//  data_set ds;
//  double maxd, maxnr, maxnb;
//  clock_t c1, c2;
//
//  if (argc > 2) print_usage();
//  if (argc == 2)  n_pts = atoi(argv[1]);
//
//  random_Gaussian_ds(&ds, n_pts);
//
//  c1 = clock();
//  set_red_blue_weight_ds(&ds, red_weight, blue_weight);
//  create_heap_ds(&ds);
//  maxd = max_discrepancy_ds(&ds, &xmin, &xmax, &ymin, &ymax, &maxnr, &maxnb);
//  c2 = clock();
//
//  printf ("maximum discrepancy: %f -- time:%f    [%1.2f,%1.2f]x[%1.2f,%1.2f]\tR:%1.2f, B:%1.2f\n",
//	  maxd,  (double)(c2-c1)/CLOCKS_PER_SEC, xmin, xmax, ymin, ymax, maxnr, maxnb);
//
//  return 0;
//}
