//
// Created by mmath on 2/12/18.
//
#include <tuple>
#include <algorithm>
#include <cstddef>
#include <limits>

#include "PolyTree.hpp"

namespace pyscan {

	Point<> const& Segment::get_line() const {
		return *this;
	}

	Point<> const& Segment::get_left() const {
	    return l_end_pt;
	}

    Point<> const& Segment::get_right() const {
        return r_end_pt;
    }


	bool Segment::lte(Point<> const& line) const {
        if (parallel(*this, line)) {
            double a = this->normalization_factor(line);
            return alte(line[2], this->operator[](2) * a);
        } else {
            return line.above_closed(l_end_pt) && line.above_closed(r_end_pt);
        }
	}

    bool Segment::lt(Point<> const& line) const {
        if (parallel(*this, line)) {
            double a = this->normalization_factor(line);
            return alt(line[2], this->operator[](2) * a);
        } else {
            return line.above(l_end_pt) && line.above(r_end_pt);
        }
    }


    bool Segment::gt(Point<> const& line) const {
        if (parallel(*this, line)) {
            double a = this->normalization_factor(line);
            return alt(this->operator[](2) * a, line[2]);
        } else {
            return line.below(l_end_pt) && line.below(r_end_pt);
        }
    }

    bool Segment::gte(Point<> const& line) const {
        if (parallel(*this, line)) {
            double a = this->normalization_factor(line);
            return alte(this->operator[](2) * a, line[2]);
        } else {
            return line.below_closed(l_end_pt) && line.below_closed(r_end_pt);
        }
    }

    bool Segment::crossed(Point<> const& line) const {
        return (line.above(l_end_pt) && line.below(r_end_pt)) ||
               (line.below(l_end_pt) && line.above(r_end_pt));
	}

	/*
	 * Returns tuple where first element is above seg and second is below seg.
	 */
    std::tuple<Segment, Segment> Segment::split(Point<> const& seg) const {
	    auto end_pt = intersection(*this, seg);
	    auto seg1 = Segment(*this, this->l_end_pt, end_pt);
	    auto seg2 = Segment(*this, end_pt, this->r_end_pt);
	    if (seg.above_closed(this->l_end_pt) && seg.below_closed(this->r_end_pt)) {
            return std::make_tuple(seg2, seg1);
	    } else {
            return std::make_tuple(seg1, seg2);
        }
	}

    /*
     * Returns tuple where first element is above seg and second is below seg.
     */
    std::tuple<WSegment, WSegment> WSegment::wsplit(Point<> const& seg) const {
        auto end_pt = intersection(*this, seg);
        if (seg.above_closed(this->get_left()) && seg.below_closed(this->get_right())) {
            return std::make_tuple(WSegment(*this, end_pt, this->get_right(), this->get_weight()),
                    WSegment(*this, this->get_left(), end_pt, this->get_weight()));
        } else {
            return std::make_tuple(WSegment(*this, this->get_left(), end_pt, this->get_weight()),
                    WSegment(*this, end_pt, this->get_right(), this->get_weight()));
        }
    }

    std::tuple<std::shared_ptr<PolyNode>, std::shared_ptr<PolyNode>>
    PolyNode::split(Point<> const& line) {
	    /*
	     * Line the left pt is the
	     */
	    Segment seg = line_to_segment(seg);
    	auto p1 = std::make_shared<PolyNode>();
    	auto p2 = std::make_shared<PolyNode>();

    	/*
    	 * Split the boundary segments
    	 */
        segment_list upper_boundary_seg;
        segment_list lower_boundary_seg;

        for (auto bs_it = boundary_seg.begin(); bs_it != boundary_seg.end(); ++bs_it) {
            if (seg.crossed(*bs_it)) {
                Segment upper, lower;
                std::tie(upper, lower) = bs_it->split(seg);
                lower_boundary_seg.push_back(lower);
                upper_boundary_seg.push_back(upper);
            } else if (seg.gte(*bs_it)) {
                upper_boundary_seg.push_back(*bs_it);
            } else {
                lower_boundary_seg.push_back(*bs_it);
            }
        }
        upper_boundary_seg.push_back(seg);
        lower_boundary_seg.push_back(seg);

        /*
         * split the crossing segments
         */
        wsegment_list upper_crossing_seg;
        wsegment_list lower_crossing_seg;

        for (auto bs_it = crossing_segments.begin(); bs_it != crossing_segments.end(); ++bs_it) {
            if (seg.crossed(*bs_it)) {
                WSegment upper, lower;
                std::tie(upper, lower) = bs_it->wsplit(seg);
                lower_crossing_seg.push_back(lower);
                upper_crossing_seg.push_back(upper);
            } else if (seg.gte(*bs_it)) {
                upper_crossing_seg.push_back(*bs_it);
            } else {
                lower_crossing_seg.push_back(*bs_it);
            }
        }

        /*
         * Split the points.
         */

        line_list upper_pts;
        line_list lower_pts;
        for (auto p_it = points.begin(); p_it != points.end(); ++p_it) {
            if (seg.above_closed(*p_it)) {
                upper_pts.push_back(*p_it);
            } else {
                lower_pts.push_back(*p_it);
            }
        }

        return make_tuple(
                std::make_shared<PolyNode>(upper_crossing_seg, upper_pts, upper_boundary_seg),
                std::make_shared<PolyNode>(lower_crossing_seg, lower_pts, lower_boundary_seg));
    }


    template <typename T>
    auto weighted_choice(T begin, T end, std::function<double(decltype(*begin) const&)> f) -> decltype(*begin) {
        std::uniform_real_distribution<double> unif(0, 1);
        std::default_random_engine re;
        double random_d = unif(re);

        double cum_weight = 0;
        double total_weight = std::accumulate(begin, end, 0, [&](double t, decltype(*begin) const& item) {
            return t + f(item);
        });
        for (; begin != end; ++begin) {
            if (cum_weight / total_weight >= random_d) {
                return *begin;
            }
        }
        return *(end - 1);
    }

    Segment PolyNode::good_line_split() const {
        auto wseg = weighted_choice(crossing_segments.begin(), crossing_segments.end(), [&](WSegment const& seg){
            return seg.get_weight();
        });

        return Segment(wseg.get_line(), wseg.get_left(), wseg.get_right());
    }

    Segment PolyNode::good_vertex_split() const {
        using int_t = decltype(crossing_segments.size());
        auto first_pt = boundary_seg.begin();
        auto second_pt = boundary_seg.begin() + (boundary_seg.end() - boundary_seg.begin()) / 2;
        return Segment(intersection(first_pt, second_pt), first_pt, second_pt);
    }

    double PolyNode::get_weight() const {
        return total_weight;
    }

    Segment PolyNode::line_to_segment(Point<> const& line) const {
        /*
         * This line crosses the polynode in one or two spots.
         */
        Point<> p[2];
        int i = 0;
        for (auto& seg : boundary_seg) {
            if (i >=  2) {
                break;
            }
            if (seg.crossed(line)) {
                p[i] = intersection(seg, line);
                i++;
            }
        }
        assert(i != 0);
        if (i == 2) {
            return Segment(line, p[0], p[1]);
        } else {
            return Segment(line, p[0], p[0]);
        }
    }


    std::tuple<Segment, Segment> four_point_split(double eps) const {
        //TODO
        return;
    }

    node_list cutting_mc(int b, node_ptr p) {

        auto by_weight = [](node_ptr const& p1, node_ptr const& p2) {
            return p1->get_weight() > p2.get_weight();
        };
        std::priority_queue<node_ptr, std::vector<node_ptr>, decltype(by_weight)> curr_nodes(by_weight);
        curr_nodes.push(p);

        while (curr_nodes.size() <= b) {
            auto& tp_el = curr_nodes.top();
            if (tp_el->get_weight() <= 0) {
                //TODO fix this line
                return node_list(curr_nodes);
            }
            //TODO finish this.
        }
    };
};