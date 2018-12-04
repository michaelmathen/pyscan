//
// Created by mmath on 2/12/18.
//
#include <tuple>
#include <algorithm>
#include <cstddef>
#include <limits>

#include "Segment.hpp"

namespace pyscan {

	Point<> Segment::get_e1() const {
	    return l_end_pt;
	}

    Point<> Segment::get_e2() const {
        return r_end_pt;
    }


	bool Segment::lte(Point<> const& loc_line) const {
        if (is_parallel(line, loc_line)) {
            return line.parallel_lte(loc_line);
        } else {
            return loc_line.above_closed(l_end_pt) && loc_line.above_closed(r_end_pt);
        }
	}

    bool Segment::lt(Point<> const& loc_line) const {
        if (is_parallel(line, loc_line)) {
            return line.parallel_lt(loc_line);
        } else {
            return loc_line.above(l_end_pt) && loc_line.above(r_end_pt);
        }
    }


    bool Segment::gt(Point<> const& loc_line) const {
        if (is_parallel(line, loc_line)) {
            return loc_line.parallel_lt(line);
        } else {
            return loc_line.below(l_end_pt) && loc_line.below(r_end_pt);
        }
    }

    bool Segment::gte(Point<> const& loc_line) const {
        if (is_parallel(line, loc_line)) {
            return line.parallel_lte(loc_line);
        } else {
            return loc_line.below_closed(l_end_pt) && loc_line.below_closed(r_end_pt);
        }
    }

    bool Segment::crossed(Point<> const& loc_line) const {
	    return loc_line.crosses(l_end_pt, r_end_pt);
	}

	/*
	 * Returns tuple where first element is above seg and second is below seg.
	 */
    std::tuple<Segment, Segment> Segment::split(Point<> const& seg) const {
	    auto end_pt = intersection(line, seg);
	    auto seg1 = Segment(*this, this->l_end_pt, end_pt);
	    auto seg2 = Segment(*this, end_pt, this->r_end_pt);
	    //std::cout << seg1 << " " << std::endl;
	    if (seg.above_closed(this->l_end_pt) && seg.below_closed(this->r_end_pt)) {
            return std::make_tuple(seg2, seg1);
	    } else {
            return std::make_tuple(seg1, seg2);
        }
	}

};