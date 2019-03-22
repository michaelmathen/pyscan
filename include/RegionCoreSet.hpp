//
// Created by mmath on 2/28/19.
//

#ifndef PYSCAN_REGIONCORESET_HPP
#define PYSCAN_REGIONCORESET_HPP

#include "Point.hpp"
#include "TrajectoryCoreSet.hpp"


namespace pyscan {


    /*
     * This takes a point list that represents a polygon and then grids the polygon.
     */
    point_list_t polygon_grid(point_list_t const& pts, double grid_r);

    point_list_t polygon_grid_hull(point_list_t pts, double min_radius, double alpha);

    point_list_t polygon_grid_even(point_list_t pts, double grid_resolution, double boundary_resolution);


    wpoint_list_t polygon_sample(const std::vector<point_list_t>& polygons, const std::vector<double>& weights,
                                 size_t sample_size);
}


#endif //PYSCAN_REGIONCORESET_HPP
