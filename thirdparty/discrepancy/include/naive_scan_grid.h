//
// Created by mmath on 3/24/19.
//

#ifndef PYSCAN_NAIVE_SCAN_GRID_H
#define PYSCAN_NAIVE_SCAN_GRID_H
#include "dataset.h"
#include "time.h"
#include "function.h"

double naive_find_rect_grid (data_set *ds, d_fxn *f,
                             int grid_size, double **red, double **blue,
                             double *xmin, double *xmax, double *ymin, double *ymax);

#endif //PYSCAN_NAIVE_SCAN_GRID_H
