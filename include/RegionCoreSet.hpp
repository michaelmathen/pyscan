//
// Created by mmath on 2/28/19.
//

#ifndef PYSCAN_REGIONCORESET_HPP
#define PYSCAN_REGIONCORESET_HPP
#include <boost/geometry/geometries/adapted/boost_polygon/polygon.hpp>

#include "Point.hpp"


namespace pyscan {

    class Polygon {
    public:
        Polygon();
        Polygon(point_list_t vertices);

        //Checks to see if the point is contained inside the polygon.
        void contains(pt2_t const& pt);


    protected:
        point_list_t vertices;
    };
}

#endif //PYSCAN_REGIONCORESET_HPP
