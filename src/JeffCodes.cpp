//
// Created by mmath on 3/24/19.
//

#include "JeffCodes.hpp"

extern "C"
{
#include "dataset.h"
#include "function.h"
#include "naive_scan.h"
#include "scan_grid.h"
#include "naive_scan_grid.h"
#include "stat_rect.h"
}


namespace pyscan {

    void set_ds(data_set* ds, wpoint_list_t const &red, wpoint_list_t const &blue) {
        ds->npts = static_cast<int>(red.size() + blue.size());
        ds->data = new data_point[ds->npts];

        size_t i = 0;
        for (auto& r : red) {
            ds->data[i].x = r[0];
            ds->data[i].y = r[1];
            ds->data[i].red = r.get_weight();
            ds->data[i].blue = 0;
        }
        for (auto& b : blue) {
            ds->data[i].x = b[0];
            ds->data[i].y = b[1];
            ds->data[i].blue = b.get_weight();
            ds->data[i].red = 0;
        }
    }

    void free_ds(data_set* ds) {
        delete[] ds->data;
        free(ds->iheap);
        free(ds->xorder);
        free(ds->yorder);
    }

    void set_func(ad_fxn* A, d_fxn* D, F_Type func, double eps, int npts) {
        switch (func) {
            case F_Type::BERNOULLI:
                set_Bernoulli_fxn(D);
                break;
            case F_Type::GAMMA:
                set_gamma_fxn(D);
                break;
            case F_Type::GAUSSIAN:
                set_Gaussian_fxn(D);
                break;
            case F_Type::POISSON:
                set_Poisson_fxn(D);
                break;
        }
        smart_gridding_fxn(D, A, eps, npts);
    }

    void set_func(d_fxn* D, F_Type func) {
        switch (func) {
            case F_Type::BERNOULLI:
                set_Bernoulli_fxn(D);
                break;
            case F_Type::GAMMA:
                set_gamma_fxn(D);
                break;
            case F_Type::GAUSSIAN:
                set_Gaussian_fxn(D);
                break;
            case F_Type::POISSON:
                set_Poisson_fxn(D);
                break;
        }
    }

    void free_func(ad_fxn* A) {
        free(A->x);
        free(A->y);
        free(A->z);
    }

    std::tuple<Rectangle, double> naive_find_rect(
            wpoint_list_t const &red,
            wpoint_list_t const &blue,
            F_Type func){

        data_set ds;
        d_fxn D;
        double xmin, xmax, ymin, ymax;
        set_func(&D, func);
        set_ds(&ds, red, blue);
        create_heap_ds(&ds);
        double m = naive_find_rect(&ds, &D, &xmin, &xmax, &ymin, &ymax);
        free_ds(&ds);
        return std::make_tuple(Rectangle(xmax, ymax, xmin, ymin), m);
    }

    std::tuple<Rectangle, double> naive_scan_grid(
            int grid_size,
            wpoint_list_t const &red,
            wpoint_list_t const &blue,
            F_Type func){
        double **red_vals;
        double **blue_vals;
        double minx, miny, dx, dy;
        data_set ds;
        d_fxn D;
        double xmin, xmax, ymin, ymax;
        set_func(&D, func);
        set_ds(&ds, red, blue);
        create_heap_ds(&ds);

        load_grid_ds(&ds, grid_size, &red_vals, &blue_vals, &minx, &dx, &miny, &dy);
        double m = naive_find_rect_grid(&ds, &D, grid_size, red_vals, blue_vals, &xmin, &xmax, &ymin, &ymax);

        xmin = xmin*dx + minx;
        xmax = xmax*dx + minx;
        ymin = ymin*dy + miny;
        ymax = ymax*dy + miny;

        free_ds(&ds);
        return std::make_tuple(Rectangle(xmax, ymax, xmin, ymin), m);
    }

    std::tuple<Rectangle, double> scan_grid(
            int grid_size,
            wpoint_list_t const &red,
            wpoint_list_t const &blue,
            double eps,
            F_Type func){
        double **red_vals;
        double **blue_vals;
        double minx, miny, dx, dy;
        data_set ds;
        d_fxn D;
        ad_fxn A;
        double m, maxd=0;

        double x_min, x_max, y_min, y_max;
        double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
        set_ds(&ds, red, blue);
        set_func(&A, &D, func, eps, ds.npts);

        smart_gridding_fxn(&D, &A, eps, ds.npts);
        load_grid_ds(&ds, grid_size, &red_vals, &blue_vals, &minx, &dx, &miny, &dy);
        create_grid_heap_ds(&ds, grid_size);

        // naively check max discrepancy
        for (int i = 0; i < A.n; ++i) {
            set_red_blue_weight_ds(&ds, A.x[i]/ds.rpts, A.y[i]/ds.bpts);
            double m = max_grid_discrepancy_ds(&ds, grid_size, red_vals, blue_vals,
                                        &x_min, &x_max, &y_min, &y_max) + A.z[i];
            if (m > maxd) {
                maxd = m;
                xmin = x_min;
                xmax = x_max;
                ymin = y_min;
                ymax = y_max;
            }
        }

        xmin = xmin*dx + minx;
        xmax = xmax*dx + minx;
        ymin = ymin*dy + miny;
        ymax = ymax*dy + miny;

        free_ds(&ds);
        free_func(&A);
        return std::make_tuple(Rectangle(xmax, ymax, xmin, ymin), m);
    }

    std::tuple<Rectangle, double> naive_approx_find_rect(
            wpoint_list_t const &red,
            wpoint_list_t const &blue,
            double eps,
            F_Type func){

        data_set ds;
        ad_fxn A;
        d_fxn D;
        double xmin, xmax, ymin, ymax;
        set_ds(&ds, red, blue);
        set_func(&A, &D, func, eps, ds.npts);
        create_heap_ds(&ds);
        double m = naive_approx_find_rect(&ds, &A, &xmin, &xmax, &ymin, &ymax);
        free_ds(&ds);
        free_func(&A);
        return std::make_tuple(Rectangle(xmax, ymax, xmin, ymin), m);
    }



    std::tuple<Rectangle, double> find_rect(
            wpoint_list_t const &red,
            wpoint_list_t const &blue,
            double eps,
            F_Type func){

        data_set ds;
        ad_fxn A;
        d_fxn D;
        double xmin=0, xmax=0, ymin=0, ymax=0, maxnr, maxnb;
        double x_min, x_max, y_min, y_max;
        double m, maxd=0;

        set_ds(&ds, red, blue);
        set_func(&A, &D, func, eps, ds.npts);
        create_heap_ds(&ds);
        for (int i = 0; i < A.n; ++i) {
            set_red_blue_weight_ds(&ds, A.x[i]/ds.rpts, A.y[i]/ds.bpts);
            m = max_discrepancy_ds(&ds, &x_min, &x_max, &y_min, &y_max, &maxnr, &maxnb) + A.z[i];
            if (m > maxd) {
              maxd = m;
              xmin = x_min;
              xmax = x_max;
              ymin = y_min;
              ymax = y_max;
            }
        }
        free_ds(&ds);
        free_func(&A);
        return std::make_tuple(Rectangle(xmax, ymax, xmin, ymin), maxd);
    }
}