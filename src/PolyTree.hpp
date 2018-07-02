//
// Created by mmath on 2/12/18.
//

#ifndef PYSCAN_CUTTINGS_HPP
#define PYSCAN_CUTTINGS_HPP

#include <vector>
#include <memory>

#include "Point.hpp"

namespace pyscan {



    class Segment : Point<2> {
        Point<> l_end_pt;
        Point<> r_end_pt;
    public:
        Segment(Point<> const& line, Point<> const& el, Point<> const& er) : Point<2>(line),
            l_end_pt(el), 
            r_end_pt(er) {
        }
        Segment() {}

        Point<> get_left() const;
        Point<> get_right() const;
        Point<> get_line() const;
        bool lte(Point<> const& line) const;
        bool gte(Point<> const& line) const;
        bool lt(Point<> const& line) const;
        bool gt(Point<> const& line) const;
        bool crossed(Point<> const& line) const;
        std::tuple<Segment, Segment> split(Point<> const& line) const;
    };

    class WSegment : Segment {
        double weight;
    public:
        WSegment(WPoint<> const& line,
                 Point<> const& el,
                 Point<> const& er, double w) : Segment<>(line, el, er),
                     weight(w) {
        }
        WSegment() : Segment<>(), weight(0){}

        double get_weight() const {
            return w;
        }
    };

    class Node {

    public:
        bool is_terminal() const;
    };

    using node_ptr = std::shared_ptr<Node>;
    using line_list = std::vector<WPoint<>>;
    using wsegment_list = std::vector<WSegment>;
    using segment_list = std::vector<Segment>;
    using node_list = std::vector<node_ptr>;


    class SegmentNode : public Node {
        node_ptr upper;
        node_ptr lower;
        Segment line;
    public:
    };

    class PolyNode : public Node {
        segment_list boundary_seg;
        wsegment_list crossing_segments;
        line_list points;
        double total_weight;

        PolyNode(wsegment_list const& ls, line_list const& pts, segment_list const& bd_seg) :
                boundary_seg(bd_seg),
                crossing_segments(ls),
                points(pts)
        {
            total_weight = std::accumulate(weights.begin(), weights.end(), 0);
        }

    public:

        /*
         * Create a polynode inside a 0,1, 0,1 box
         */
        PolyNode(wsegment_list const& cs, segment_list const& bs, line_list const& pts) : points(pts),
                                                                                          crossing_segments(cs),
                                                                                          boundary_seg(bs)
        {


            boundary_seg = { {intersection({0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}), {0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}},
                             {intersection({1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}), {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}},
                             {intersection({1.0, 1.0, 1.0}, {1.0, 0.0, 1.0}), {1.0, 1.0, 1.0}, {1.0, 0.0, 1.0}},
                             {intersection({1.0, 0.0, 1.0}, {0.0, 0.0, 1.0}), {1.0, 0.0, 1.0}, {0.0, 0.0, 1.0}}};

            total_weight = std::accumulate(crossing_segments.begin(), crossing_segments.end(), 0, [&](double w, WSegment const& seg){
                return w + seg.get_weight();
            });
        }

        std::tuple<std::shared_ptr<PolyNode>, std::shared_ptr<PolyNode>>
        split(Point<> const& seg);

        Segment line_to_segment(Point<2> const& ) const;
        Segment good_line_split() const;
        Segment good_vertex_split() const;

        /*
         * Find a split that approximately splits the points into 4 segments with approximately the same points.
        */
        Segment four_point_split(double eps) const;

        /*
         * ham sandwich split till max weight
        */
        node_list partition_mw(double max_weight, double eps);

        /*
         * ham sandwich split till there are b cells.
        */
        node_list partition_mc(int b);

        /*
         *Cutting till max weight
        */
        node_list cutting_mw(double max_weight);

        /*
         * Cut till there are b cells.
        */
        node_list cutting_mc(int b);


        double get_weight() const;

    };

}
#endif //PYSCAN_CUTTINGS_HPP
