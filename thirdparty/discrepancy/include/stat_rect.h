//
// Created by mmath on 3/24/19.
//

#ifndef PYSCAN_STAT_RECT_H
#define PYSCAN_STAT_RECT_H

#include "dataset.h"
#include "function.h"

double naive_approx_find_rect (data_set *ds, ad_fxn *A,
                               double *xmin, double *xmax, double *ymin, double *ymax);

#endif //PYSCAN_STAT_RECT_H
