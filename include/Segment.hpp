//
// Created by mmath on 2/12/18.
//

#ifndef PYSCAN_SEGMENT_HPP
#define PYSCAN_SEGMENT_HPP

#include <vector>
#include <memory>

#include "Point.hpp"

namespace pyscan {



    class Segment : public Point<2> {
        Point<> l_end_pt;
        Point<> r_end_pt;
    public:

        explicit Segment(Point<2> const& line) : Point<2>(line.orient_down()),
                    l_end_pt(-line[1], line[0], 0.0),
                    r_end_pt(line[1], -line[0], 0.0) {}

        Segment(Point<2> const& line, Point<2> const& el, Point<2> const& er) : Point<2>(line.orient_down()),
            l_end_pt(el), 
            r_end_pt(er) {
        }
        Segment() {}

        Point<> get_e1() const;
        Point<> get_e2() const;
        bool lte(Point<> const& line) const;
        bool gte(Point<> const& line) const;
        bool lt(Point<> const& line) const;
        bool gt(Point<> const& line) const;
        bool crossed(Point<> const& line) const;
        std::tuple<Segment, Segment> split(Point<> const& seg) const;

        std::string str() const override {
            std::stringstream ss;
            ss << "Segment(" << *this;
            ss << ", " << this->get_e1();
            ss << ", " << this->get_e2();
            ss << ")";
            return ss.str();
        }
    };

}
#endif
