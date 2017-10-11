//
// Created by mmath on 7/26/17.
//

#ifndef PYSCAN_EPSSAMPLES_HPP
#define PYSCAN_EPSSAMPLES_HPP
#include <tuple>
#include <vector>

#include "Point.hpp"

namespace pyscan {
    enum class TrapType {
        NORMAL,
        L_Deg,
        R_Deg
    };

    class Trapezoid;

    using Line = std::tuple<double, double, int>; // x, y, crossing number
    using Line2 = std::tuple<double, double>;
    using PointList = std::vector<Line2>;
    using LineList = std::vector<Line>;
    using Cell = std::tuple<Trapezoid, LineList, PointList, long>; // Trapezoid, lines, points, weight

    void splitCell(Cell const&, Line const&, std::vector<Cell>&);

    class Trapezoid {
        // If degenerate is true than x1, y1, x2, y2, x3, y3 define a triangle.
    public:
        double top_a, top_b, bottom_a, bottom_b, left_x, right_x;
        TrapType degenerate = TrapType::NORMAL;
        Trapezoid();
        Trapezoid(double ta, double tb, double lx, double ba, double bb, double rx);
        bool crossing(double a, double b) const;
        bool incident(double a, double b) const;
        bool crossesInterior(double a, double b) const;

        void setBottomLine(double a, double b);
        void setTopLine(double a, double b);

        bool crossesTop(double a, double b) const;
        bool crossesBottom(double a, double b) const;
        bool crossesLeftSide(double a, double b) const;
        bool crossesRightSide(double a, double b) const;
        bool isDegenerate() const;
        bool lDegenerate() const;
        bool rDegenerate() const;
        bool finiteTop() const;
        bool finiteBottom() const;
        bool finiteLeft() const;
        bool finiteRight() const;
        bool contains(Line2 const& pt) const;
        Line2 ulCorner() const;
        Line2 urCorner() const;
        Line2 llCorner() const;
        Line2 lrCorner() const;


        void setLeft(double left);
        void setRight(double right);

        std::string print() const;


        friend void splitCell(Cell const&, Line const&, std::vector<Cell>&);

    };


    LineList dualizeToLines(std::vector<Cell> const& cells);
    std::vector<Cell> createPartition(PointList & points, int r);

    // Defines functions for creating small size eps-samples for
    // rectangles and halfpspaces.
    using point_it= std::vector<Point<>>::iterator;

    void rectangleSample(point_it begin, point_it end,
                         point_it output_b, point_it output_e);

    /*
     * void halfpsaceSample(point_it begin, point_it end,
     *                   point_it output_b, point_it output_e);
     */
    std::vector<Cell> randomCuttings(LineList &, PointList &, int);

    double xLoc(double, double, double, double);
    double xLoc(Line2, Line2);

    double yVal(Line2, double);

//    class Edge;
//    class Simplex;
//    class Vertex;
//
//    using edge_ptr = std::shared_ptr<Edge>;
//    using simplex_ptr = std::shared_ptr<Simplex>;
//    using vertex_ptr = std::shared_ptr<Vertex>;
//
//    class Edge {
//    public:
//        Line2 edge_line;
//        double x_begin;
//        double x_end;
//        std::weak_ptr<Simplex> bottom_ptr;
//        std::weak_ptr<Simplex> top_ptr;
//
//        Edge(Line2 l, double x1, double x2, simplex_ptr p_top, simplex_ptr p_bot) :
//                edge_line(std::move(l)),
//                x_begin(std::min(x1, x2)),
//                x_end(std::max(x1, x2)),
//                bottom_ptr(p_bot),
//                top_ptr(p_top)
//                {}
//
//        Edge(Edge old_edge, double x1, double x2) :
//                Edge(old_edge.edge_line, x1, x2, old_edge.top_ptr.lock(), old_edge.bottom_ptr.lock()) {
//        }
//        bool crosses(Line2 const& line) const;
//        bool sameLine(Edge const& other_edge) const;
//
//        bool posSlope() const {
//            return std::get<0>(edge_line) > 0;
//        }
//    };
//
//    class Vertex {
//        /*
//         * Weirdly two vertices can be on top of each other.
//         * This just provides a mapping where two lines cross each other.
//         */
//    public:
//        std::weak_ptr<Simplex> left_ptr;
//        std::weak_ptr<Simplex> right_ptr;
//        Line2 vertex_location;
//        Vertex(Line2 const& l1, Line2 const& l2, simplex_ptr pl, simplex_ptr pr);
//        bool crosses(Line2 const& line) const;
//    };
//
//
//
//    class Splitter {
//        /*
//         * The point where line_left and line_right meet define
//         * a vertical line. The line oppo defines the height of the splitter
//         *  line_left, line_right_
//         *  _
//         * /|
//         *  |
//         *  |/
//         *  /
//         * line_oppo
//         *
//         * Need to be able to quickly tell if these splitters are still valid.
//         */
//        bool _standing;
//        bool _unbounded = false;
//        void setLines(Line2 const& l1, Line2 const& l2) {
//            if (_standing) {
//                if (std::get<0>(l1) <= std::get<0>(l2)) {
//                    l = l1, r = l2;
//                } else {
//                    l = l2, r = l1;
//                }
//            } else {
//                if (std::get<0>(l1) <= std::get<0>(l2)) {
//                    l = l2, r = l1;
//                } else {
//                   l = l1, r = l2;
//                }
//            }
//        }
//    public:
//        Line2 l;
//        Line2 r;
//        Line2 termination;
//
//        Splitter(Line2 const& l1, Line2 const& l2, bool standing) : _standing(standing), _unbounded(true) {
//            setLines(l1, l2);
//        }
//
//        Splitter(Line2 l1, Line2 l2, Line2 term)  {
//            double x_val = getX();
//            _standing = yVal(l, x_val) <= yVal(term, x_val);
//            setLines(l1, l2);
//        }
//        Splitter(Splitter old_splitter, Line2 t) : Splitter(old_splitter), termination(t) {}
//
//        double getX() const;
//        double bottom() const;
//        double top() const;
//
//        bool crosses(Line2 const& l) const;
//
//        bool crosses(Line2 const& l, Line2& pt) const;
//
//        bool crosses(Line const& l) const;
//
//        bool standing() const { return _standing; }
//        bool passesBelow(Line2 const& l) const;
//        bool passesAbove(Line2 const& l) const;
//        bool leftOf(Line2 const& l) const;
//        bool rightOf(Line2 const& l) const;
//        Line2 crossPt() const;
//        bool isCrossLine(Line2 const& nl) const;
//        bool leftOf(Line const& l) const;
//        bool rightOf(Line const& l) const;
//    };
//
//    void splitSimplex(simplex_ptr& sp, simplex_ptr& entr_ptr, simplex_ptr& upper, simplex_ptr& lower);
//
//    class Simplex {
//        /*
//         * Contains
//         */
//        bool open = true;
//        simplex_ptr this_ptr;
//        // Defines the edges of the simplex.
//        //These are stored in counter clockwise order
//        // If the simplex is open then the opening lies between the begining and ending of the edge list.
//        std::vector<edge_ptr> edge_list;
//        // Pointers to the simplices that are adjacent to the corners of this simplex.
//        // The ith pointer lies at the beginning of the ith edge.
//        std::vector<vertex_ptr> corner_ptrs;
//
//        // These define vertical lines that subdivide the inside
//        // of the simplex into a vertical decomposition.
//        // They are stored in x increasing order.
//        std::vector<Splitter> splitter_list;
//        // The list of lines that cross each region defined by the splitters.
//        // There will always be n+1 of these where n is the number of splitters.
//        // They are stored in x increasing order.
//        std::vector<LineList> crossings;
//
//        /*
//         * Splits the splitters based on the new splitting line.
//         * Does not merge splitters that are now defined by line_left
//         * and line_right that are outside of this simplex.
//         */
//
//        Splitter makeSplitter(edge_ptr e1, Line2 new_line, edge_ptr term);
//        void splitSplitters(Line2 const& new_line, simplex_ptr entr_ptr,
//                            decltype(crossings)& above_crossings, decltype(splitter_list)& above_splitters,
//                            decltype(crossings)& below_crossings, decltype(splitter_list)& below_splitters);
//
//
//        void copyTopHalf(int entrance, int exit,
//                         edge_ptr entrance_edge,
//                         edge_ptr exit_edge,
//                         edge_ptr cutting_edge,
//                         vertex_ptr entrance_v,
//                         vertex_ptr exit_v,
//                         std::vector<edge_ptr>& top_edges,
//                         std::vector<vertex_ptr>& top_vertices);
//
//        void copyBottomHalf(int entrance, int exit,
//                            edge_ptr entrance_edge,
//                            edge_ptr exit_edge,
//                            edge_ptr cutting_edge,
//                            vertex_ptr entrance_v,
//                            vertex_ptr exit_v,
//                            std::vector<edge_ptr>& top_edges,
//                            std::vector<vertex_ptr>& top_vertices);
//
//        edge_ptr findOppositeEdge(edge_ptr intersecting_edge, double x_val);
//        void computeLRVertexCross(simplex_ptr left_ptr, Line2 l, vertex_ptr left_vt, vertex_ptr right_vt);
//        void computeLRcrossing(simplex_ptr left_ptr, Line2 l, edge_ptr left_edge, edge_ptr right_edge) const;
//        void updatePtrs(simplex_ptr);
//    public:
//        Simplex upper_split(Line2 const& new_line, edge_ptr clock_wise, vertex_ptr lower_simplex);
//        bool closed() const;
//    };
//
//    class Arrangement {
//    public:
//
//    };

}


#endif //PYSCAN_EPSSAMPLES_HPP
