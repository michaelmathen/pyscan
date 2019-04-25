//
// Created by mmath on 2/28/19.
//
#include <numeric>

#include "RegionCoreSet.hpp"
#include "Sampling.hpp"
#include "TrajectoryCoreSet.hpp"

#include <functional>


#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Triangle_2.h>

namespace pyscan {

    struct FaceInfo2
    {
        FaceInfo2(){}
        int nesting_level;
        bool in_domain(){
            return nesting_level%2 == 1;
        }
    };

    using K                      = CGAL::Simple_cartesian<double>;
    using Traits                 = CGAL::Partition_traits_2<K>;
    using cgal_pt_t              = Traits::Point_2;
    using cgal_poly_t            = Traits::Polygon_2;
    using vertex_iterator_t      = cgal_poly_t::Vertex_iterator;
    using cgal_poly_list_t       = std::vector<cgal_poly_t>;
    using Creator                = CGAL::Creator_uniform_2<int, cgal_pt_t>;
    using Point_generator        = CGAL::Random_points_in_square_2<cgal_pt_t, Creator >;


    typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K>    Fbb;
    typedef CGAL::Exact_predicates_tag                                Itag;
    typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb> Fb;
    typedef CGAL::Triangulation_vertex_base_2<K> Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;


    inline static std::tuple<double, double, double, double> bounding_box(point_list_t::const_iterator traj_b,
                                                            point_list_t::const_iterator traj_e) {

        double lx = (*traj_b)(0), ux = (*traj_b)(0), ly = (*traj_b)(1), uy = (*traj_b)(1);
        for ( ; traj_b != traj_e; traj_b++) {
            lx = std::min(lx, (*traj_b)(0));
            ux = std::max(ux, (*traj_b)(0));
            ly = std::min(ly, (*traj_b)(1));
            uy = std::max(uy, (*traj_b)(1));
        }
        return std::make_tuple(lx, ly, ux, uy);
    }


    inline static double y_to_x(double x_1, double y_1, double x_2, double y_2, double y) {
        return (x_1 - x_2) / (y_1 - y_2) * (y - y_1) + x_1;

    }


    inline static bool inside_box(double x, double y, pt2_t const& lpt, pt2_t const& rpt) {
        double ux = std::max(lpt(0), rpt(0));
        double lx = std::min(lpt(0), rpt(0));
        double uy = std::max(lpt(1), rpt(1));
        double ly = std::min(lpt(1), rpt(1));
        return util::alte(x, ux) && util::alte(lx, x) &&  util::alte(y, uy) && util::alte(ly, y);
    }

    std::vector<std::vector<Point<2>>> grid_horz(point_list_t::const_iterator traj_b,
                                                 point_list_t::const_iterator traj_e,
                                                 double chord_l,
                                                 double& ux, double &uy, double& lx, double& ly) {
        /*
         * Creates a grid where every time we cross a y division we place a point.
         */


        if (traj_e - traj_b <= 1) {
            return {};
        }

        auto last_pt = traj_e - 1;


        std::tie(lx, ly, ux, uy) = bounding_box(traj_b, traj_e);

        long g_size = static_cast<long>((uy - ly) / chord_l) + 1;
        std::vector<std::vector<Point<>>> traj_points(g_size, std::vector<pt2_t>());

        for (auto curr_pt = traj_b; curr_pt != traj_e; curr_pt++) {

            auto g_x = static_cast<int>(((*last_pt)(0) - lx) / chord_l);
            auto g_y = static_cast<int>(((*last_pt)(1) - ly) / chord_l);
            auto g_n_x = static_cast<int>(((*curr_pt)(0) - lx) / chord_l);
            auto g_n_y = static_cast<int>(((*curr_pt)(1) - ly) / chord_l);

            if (g_n_x < g_x) {
                std::swap(g_n_x, g_x);
            }
            if (g_n_y < g_y) {
                std::swap(g_n_y, g_y);
            }
            for (int j = g_y; j <= g_n_y; j++) {
                double y_val = j * chord_l + ly;
                double x_val = y_to_x((*last_pt)(0), (*last_pt)(1), (*curr_pt)(0), (*curr_pt)(1), y_val);
                if (!inside_box(x_val, y_val, *last_pt, *curr_pt)) {
                    continue;
                }
                traj_points[j].emplace_back(x_val, y_val, 1.0);
            }
            last_pt = curr_pt;
        }
        return traj_points;
    }
    /*
     * This takes a point list that represents a polygon and then grids the polygon.
     */
    point_list_t polygon_grid(point_list_t const& pts, double grid_r) {

        double ux, uy, lx, ly;

        std::vector<pt2_t> grid_pts;
        // This inserts a point wherever this crosses row.
        auto rows = grid_horz(pts.begin(), pts.end(), grid_r, ux, uy, lx, ly);

        for (auto& row : rows) {
            std::sort(row.begin(), row.end(), [](const pt2_t& p1, const pt2_t& p2) {
               return p1(0) < p2(0);
            });
            // Now that it is sorted grid the horz
            auto lpt = row.begin();
            bool inside = true;
            for (auto cpt = row.begin() + 1; cpt != row.end(); cpt++) {
                if (inside) {

                    auto g_x = static_cast<int>(((*lpt)(0) - lx) / grid_r);
                    auto g_n_x = static_cast<int>(((*cpt)(0) - lx) / grid_r);
                    for (int j = g_x + 1; j <= g_n_x; j++) {
                        double x_val = j * grid_r + lx;
                        grid_pts.emplace_back(x_val, (*lpt)(1), 1.0);
                    }
                }
                lpt = cpt;
                inside = !inside;
            }
        }
        remove_duplicates(grid_pts);

        return grid_pts;
    }



    point_list_t polygon_grid_hull(point_list_t pts, double min_radius, double alpha) {
        /*diag = 2 * r so (2 * r)^2 = 2 * gr^2 => 2 * r^2 = gr^2 => sqrt(2) * r = gr*/
        if (pts.size() <= 1) {
            return pts;
        }

        point_list_t internal_pts = polygon_grid(pts, min_radius * std::sqrt(2));


        if (alpha > min_radius * std::sqrt(2)) {
            if (internal_pts.empty()) {
                internal_pts.emplace_back(pts[0]);
            }
            return internal_pts;
        } else {
            double chord_l = std::sqrt(4 * alpha * min_radius - 2 * alpha * alpha);
            pts.emplace_back(pts.front());
            point_list_t ext_pts = approx_traj_kernel_grid(pts, chord_l, alpha);
            std::copy(ext_pts.begin(), ext_pts.end(), std::back_inserter(internal_pts));
        }
        remove_duplicates(internal_pts);
        return internal_pts;
    }



    point_list_t polygon_grid_even(point_list_t pts, double grid_resolution, double boundary_resolution) {
        /*diag = 2 * r so (2 * r)^2 = 2 * gr^2 => 2 * r^2 = gr^2 => sqrt(2) * r = gr*/
        if (pts.size() <= 2) {
            return even_sample_error(pts, boundary_resolution, false);
        }
        auto internal_pts = polygon_grid(pts, grid_resolution);
        pts.emplace_back(pts.front());
        auto ext_pts = even_sample_error(pts, boundary_resolution, false);
        std::copy(ext_pts.begin(), ext_pts.end(), std::back_inserter(internal_pts));
        remove_duplicates(internal_pts);
        return internal_pts;
    }

    cgal_poly_list_t triangulate_convex(vertex_iterator_t b, vertex_iterator_t e) {

        auto init_b = b;
        b++;
        if (b == e) { return {}; }
        auto last_b = b;
        b++;
        if (b == e) { return {}; }
        cgal_poly_list_t triangle_list;
        for (; b != e; b++) {
            triangle_list.emplace_back(cgal_poly_t());
            triangle_list.back().push_back(*init_b);
            triangle_list.back().push_back(*last_b);
            triangle_list.back().push_back(*b);
            last_b = b;
        }
        return triangle_list;
    }

    template<typename V, typename URRNG>
    wpt2_t triangle_sample(V const& triangle, URRNG&& gen) {
        std::uniform_real_distribution<double> dist(0, 1);

        double r1 = dist(gen);
        double r2 = dist(gen);

        auto v0 = triangle[0];
        auto v1 = triangle[1];
        auto v2 = triangle[2];
        double x0 = v0.x(), x1 = v1.x(), x2 = v2.x();
        double y0 = v0.y(), y1 = v1.y(), y2 = v2.y();
        double a = sqrt(r1);
        return wpt2_t(1.0,
                (1 - a) * x1 + (a * (1 - r2)) * x2 + (a * r2) * x0,
                (1 - a) * y1 + (a * (1 - r2)) * y2 + (a * r2) * y0, 1.0);
    }

    void mark_domains(CDT& ct,
                 CDT::Face_handle start,
                 int index,
                 std::list<CDT::Edge>& border ) {
        if(start->info().nesting_level != -1){
            return;
        }
        std::list<CDT::Face_handle> queue;
        queue.push_back(start);
        while(! queue.empty()){
            CDT::Face_handle fh = queue.front();
            queue.pop_front();
            if(fh->info().nesting_level == -1){
                fh->info().nesting_level = index;
                for(int i = 0; i < 3; i++){
                    CDT::Edge e(fh,i);
                    CDT::Face_handle n = fh->neighbor(i);
                    if(n->info().nesting_level == -1){
                        if(ct.is_constrained(e)) border.push_back(e);
                        else queue.push_back(n);
                    }
                }
            }
        }
    }


    void mark_domains(CDT& cdt) {
        for(CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
            it->info().nesting_level = -1;
        }
        std::list<CDT::Edge> border;
        mark_domains(cdt, cdt.infinite_face(), 0, border);
        while(! border.empty()){
            CDT::Edge e = border.front();
            border.pop_front();
            CDT::Face_handle n = e.first->neighbor(e.second);
            if(n->info().nesting_level == -1){
                mark_domains(cdt, n, e.first->info().nesting_level+1, border);
            }
        }
    }

    wpoint_list_t polygon_sample(const std::vector<point_list_t>& polygons, const std::vector<double>& weights,
                                size_t sample_size) {
        if (weights.size() != polygons.size()) {
            throw "Weights and polygons must be the same length.";
        }
        cgal_poly_list_t cpolys(polygons.size(), cgal_poly_t());

        auto it_cgal = cpolys.begin();
        auto it_pyscan = polygons.begin();

        /*
         * Construct all the polygons
         */
        for (; it_cgal != cpolys.end(); it_cgal++, it_pyscan++) {
            for (auto & pt : *it_pyscan) {
                it_cgal->push_back(cgal_pt_t(pt(0), pt(1)));
            }
        }
        /*
         * Now take a sample of the weights and then we can sample from the polygons.
         */
        std::random_device rd;
        std::minstd_rand gen(rd());
        std::vector<size_t> pt_counts(weights.size());
        {
            std::vector<size_t> indices(weights.size());
            std::vector<decltype(indices.begin())> sample_indices(sample_size);

            std::iota(indices.begin(), indices.end(), 0);
            std::function<double(decltype(indices.begin()))> wf = [&](decltype(indices.begin()) i) {
                return weights[*i];
            };
            random_sample_wr(indices.begin(), indices.end(), gen, wf, sample_size, sample_indices.begin());
            /*
             * Count how many samples are needed for each polygon.
             */
            for (auto i : sample_indices) {
                pt_counts[*i] += 1;
            }
        }

        std::vector<wpt2_t> sample_pts;
        for (size_t i = 0; i < pt_counts.size(); i++) {
            if (pt_counts[i] > 0) {
                cgal_poly_list_t partition_polys;
                CDT cdt;
                cdt.insert_constraint(cpolys[i].vertices_begin(), cpolys[i].vertices_end(), true);

                mark_domains(cdt);

                std::vector<CDT::Finite_faces_iterator> sample_face;
                std::function<double(CDT::Finite_faces_iterator)> wf = [&](CDT::Finite_faces_iterator face) {
                    if (face->info().in_domain()) {
                        return std::abs(cdt.triangle(face).area());
                    }  else {
                        return 0.0;
                    }
                };
                random_sample_wr(cdt.finite_faces_begin(), cdt.finite_faces_end(), gen, wf , pt_counts[i], std::back_inserter(sample_face));

                for (auto& face : sample_face) {
                    auto pt = triangle_sample(cdt.triangle(face), gen);
                    sample_pts.emplace_back(pt);
                }
            }
        }
        return sample_pts;
    }
}