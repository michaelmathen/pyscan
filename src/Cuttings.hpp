//
// Created by mmath on 2/12/18.
//

#ifndef PYSCAN_CUTTINGS_HPP
#define PYSCAN_CUTTINGS_HPP

#include <vector>
#include <memory>

namespace pyscan {

    class ProjObject {
        double a, b, c;

    public:
        ProjObject(double a, double b, double c) : a(a), b(b), c(c) {}

        ProjObject intersection(ProjObject const& other_line) const;
        bool dual_is_below_closed(ProjObject const& pt) const;
        double evaluate(ProjObject const& pt) const;
        bool parallel(ProjObject const& line) const;
    };


    class ProjSegment : public ProjObject {
        ProjObject l_end_pt;
        ProjObject r_end_pt;
    public:
        ProjSegment(ProjObject const& line, ProjObject const& el, ProjObject const& er) : ProjObject(line), 
            l_end_pt(el), 
            r_end_pt(er) {}

        ProjObject get_line() const;

        bool above_closed() const;
        bool below_closed() const;
        bool is_crossed() const;
    };

    using node_ptr = std::shared_ptr<Node>;
    using line_list = std::vector<ProjObject>;
    using segment_list = std::vector<ProjSegment>;
    using weight_list = std::vector<double>;

    class Node {


    };

    class SegmentNode : public Node {
        node_ptr upper;
        node_ptr lower;
        ProjObject line;
    public:
        bool crossing(ProjSegment const& segment) const;
        bool above_closed(ProjSegment const& segment) const;
        bool below_closed(ProjSegment const& segment) const;

        bool above_closed_pt(ProjObject const& pt) const;
        bool below_closed_pt(ProjObject const& pt) const;
    };

    class PolyNode : public Node {
        line_list boundary_lines;
        line_list crossing_lines;
        line_list points;
        weight_list weights;
        double total_weight;
    public:
        PolyNode(line_list const& ls, weight_list const& ws, line_list const& pts) : 
            crossing_lines(ls), 
            weights(ws),
            points(pts)
        {
            total_weight = std::accumulate(weights.begin(), weights.end(), 0);
        }

        void split(ProjObject const& line, node_ptr& p1, node_ptr& p2) const;

        ProjSegment good_line_split(int k) const;
        ProjSegment good_vertex_split(int k) const;
        double get_weight() const;

    };


    std::vector<Cell> triangle_cutting(std::vector<ProjObject> const& lines, int r);


}
#endif //PYSCAN_CUTTINGS_HPP
