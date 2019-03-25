//
// Created by mmath on 3/24/19.
//

#ifndef PYSCAN_NAIVE_SCAN_H
#define PYSCAN_NAIVE_SCAN_H
#include "dataset.h"
#include "function.h"

double naive_find_rect (data_set *ds, d_fxn *f,
                        double *xmin, double *xmax, double *ymin, double *ymax);

#endif //PYSCAN_NAIVE_SCAN_H
