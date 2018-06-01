
from seidel_tree import *
import poly_tree
import unittest
import random
import matplotlib.pyplot as plt

testing_set = [Line(**{   'a': -0.866656637164839,
    'b': 0.06593941835289452}), Line(**{   'a': 0.19404522871026875,
    'b': 0.7886667631517558}), Line(**{   'a': -0.17893431540696847,
    'b': 0.3794067891374904}), Line(**{   'a': -0.8917532021685302,
    'b': 0.054217301335525825}), Line(**{   'a': 0.715824412774232,
    'b': 0.8733584769207572}), Line(**{   'a': -0.41577400211479454,
    'b': 0.13281573615549946}), Line(**{   'a': -0.679547315053992,
    'b': 0.009189296271796654}), Line(**{   'a': -0.8270605066537602,
    'b': 0.12228767149031405}), Line(**{   'a': -0.38271037022509313,
    'b': 0.44398449063278567}), Line(**{   'a': -0.4234965292114994,
    'b': 0.05413516846246624}), Line(**{   'a': 0.3803219223045814,
    'b': 0.5695273651378064}), Line(**{   'a': -0.14465571198994198,
    'b': 0.5307055816821469}), Line(**{   'a': -0.5977082177654891,
    'b': 0.13647917590614056}), Line(**{   'a': 0.39079473043470747,
    'b': 0.42917744896680987}), Line(**{   'a': -0.05758577149225108,
    'b': 0.15678313464074756}), Line(**{   'a': -0.005265034538097568,
    'b': 0.02792193676442911}), Line(**{   'a': 0.3802621360289349,
    'b': 0.9937717503230108}), Line(**{   'a': 0.18696008854001067,
    'b': 0.9111863037351035}), Line(**{   'a': -0.23308973383053566,
    'b': 0.17540714132879687}), Line(**{   'a': 0.7471973132647067,
    'b': 0.9205304538399265}), Line(**{   'a': 0.44391857810165314,
    'b': 0.5659948185840681}), Line(**{   'a': 0.6069243279587145,
    'b': 0.6328252390501302}), Line(**{   'a': -0.11250885861109172,
    'b': 0.42460092601393307}), Line(**{   'a': -0.3077652354300694,
    'b': 0.11291003644751385}), Line(**{   'a': -0.44554145868159944,
    'b': 0.03867145102222724}), Line(**{   'a': 0.47816751595135476,
    'b': 0.8734283397600994}), Line(**{   'a': -0.45268326250319735,
    'b': 0.45267349131486934}), Line(**{   'a': 0.3408278962587603,
    'b': 0.42847021072572866}), Line(**{   'a': -0.8703154308959706,
    'b': 0.10507644279585704}), Line(**{   'a': 0.4942443868247973,
    'b': 0.6619908824271351}), Line(**{   'a': -0.5162920262469274,
    'b': 0.38560284209625273}), Line(**{   'a': 0.26271853685173274,
    'b': 0.39364428859226996}), Line(**{   'a': 0.061947278573956654,
    'b': 0.5067177339522749}), Line(**{   'a': -0.3981199910873411,
    'b': 0.152886778997951}), Line(**{   'a': 0.7615896455228987,
    'b': 0.8373099315983583}), Line(**{   'a': 0.13717420290171323,
    'b': 0.14881940828695395}), Line(**{   'a': 0.08781140217528582,
    'b': 0.9372449661769702}), Line(**{   'a': 0.31129370365977316,
    'b': 0.8204313011511879}), Line(**{   'a': 0.16401308672500214,
    'b': 0.5608815778417393}), Line(**{   'a': -0.2947980355207416,
    'b': 0.5303919040310898}), Line(**{   'a': 0.7384251351572098,
    'b': 0.9489296162716477}), Line(**{   'a': -0.7773424710420901,
    'b': 0.13923312940234045}), Line(**{   'a': 0.1966622666421921,
    'b': 0.6739223688610961}), Line(**{   'a': 0.3724358679020968,
    'b': 0.8006953947425834}), Line(**{   'a': -0.18250118322396247,
    'b': 0.4968055164181834}), Line(**{   'a': 0.3109127306848579,
    'b': 0.9735931914730755}), Line(**{   'a': 0.7020528188182354,
    'b': 0.9411558034300355}), Line(**{   'a': -0.015840141941344088,
    'b': 0.664233453811992}), Line(**{   'a': 0.3127729545556547,
    'b': 0.8009259473422059}), Line(**{   'a': 0.5801333718113135,
    'b': 0.732702829324671}), Line(**{   'a': 0.4859676399691063,
    'b': 0.9417884626987422}), Line(**{   'a': 0.3548517715571854,
    'b': 0.47390878645354084}), Line(**{   'a': 0.03619698223496404,
    'b': 0.38524081823957057}), Line(**{   'a': -0.22594527428947753,
    'b': 0.19297248446602733}), Line(**{   'a': -0.012995502317622476,
    'b': 0.15703037304763}), Line(**{   'a': -0.6252498382231659,
    'b': 0.1903702182457886}), Line(**{   'a': 0.7072991472229251,
    'b': 0.8392662781390835}), Line(**{   'a': 0.1033813452713842,
    'b': 0.8589399284628827}), Line(**{   'a': 0.049643691127804646,
    'b': 0.8254201036311279}), Line(**{   'a': -0.05817087467849169,
    'b': 0.608376813919617}), Line(**{   'a': 0.1404928726114879,
    'b': 0.7083937958141401}), Line(**{   'a': -0.04850021206178401,
    'b': 0.22776234360949832}), Line(**{   'a': -0.22832469169259706,
    'b': 0.1888990885863897}), Line(**{   'a': 0.645588768419232,
    'b': 0.8381923552831446}), Line(**{   'a': -0.4378856181271715,
    'b': 0.5208263685809627}), Line(**{   'a': -0.17579301388377644,
    'b': 0.44828199614194164}), Line(**{   'a': 0.08056615185540039,
    'b': 0.887659623379861}), Line(**{   'a': 0.2666696638753986,
    'b': 0.3344353245839955}), Line(**{   'a': -0.48821204910743976,
    'b': 0.18017753918807722}), Line(**{   'a': -0.38184269113088165,
    'b': 0.12268253311114685}), Line(**{   'a': 0.5371371947262155,
    'b': 0.9255840854344357}), Line(**{   'a': -0.03266876020231846,
    'b': 0.666423767011373}), Line(**{   'a': 0.10856512209265978,
    'b': 0.5341362506539464}), Line(**{   'a': -0.3605424753944233,
    'b': 0.36722748653979154}), Line(**{   'a': 0.022092292360553545,
    'b': 0.07304888995290681}), Line(**{   'a': 0.12006897138371186,
    'b': 0.3167113824332932}), Line(**{   'a': 0.08921234641744591,
    'b': 0.6677764319032927}), Line(**{   'a': 0.04353645811117224,
    'b': 0.26565881253202195}), Line(**{   'a': 0.4865764712065329,
    'b': 0.6203017833085297}), Line(**{   'a': 0.2833558932979773,
    'b': 0.43051108572015817}), Line(**{   'a': 0.34765349546646773,
    'b': 0.9531028494478859}), Line(**{   'a': 0.7535343906097931,
    'b': 0.9452764368453376}), Line(**{   'a': 0.6055531124626453,
    'b': 0.9718205821118342}), Line(**{   'a': -0.18248074314965634,
    'b': 0.12150851252983919}), Line(**{   'a': 0.30505598351667607,
    'b': 0.5455210504889672}), Line(**{   'a': -0.5390395830541774,
    'b': 0.027083277035031106}), Line(**{   'a': -0.4122843851600505,
    'b': 0.08079269527782074}), Line(**{   'a': 0.5278021327149751,
    'b': 0.943699209623703}), Line(**{   'a': 0.49539589071058987,
    'b': 0.8891975833691151}), Line(**{   'a': -0.744002372254413,
    'b': 0.23906732041247014}), Line(**{   'a': 0.4308809186519963,
    'b': 0.9993369810076742}), Line(**{   'a': 0.42310682742420147,
    'b': 0.6593099528631732}), Line(**{   'a': 0.10546723586132944,
    'b': 0.6180696255870675}), Line(**{   'a': -0.38779494113800705,
    'b': 0.3572132320410676}), Line(**{   'a': 0.4376003584430792,
    'b': 0.7412891126716434}), Line(**{   'a': -0.852572999061909,
    'b': 0.06797998304249575}), Line(**{   'a': -0.09252787021070796,
    'b': 0.21554809832985333}), Line(**{   'a': -0.28613331031248534,
    'b': 0.20033307511238063}), Line(**{   'a': 0.06101258925163622,
    'b': 0.2139179276696056}), Line(**{   'a': 0.6702513310091363,
    'b': 0.7492456602119376})]

class LineTestingMethod(unittest.TestCase):

    def test_above_interval(self):
        l1 = Line(**{'a': -0.28002249578781524,
                     'b': 0.04670219515585938})
        l2 = Line(**{'a': -0.08642455382757375,
                     'b': 0.29438253103340406})

        x_intercept = l1.x_intercept(l2)

        self.assertFalse(l1.above_interval(l2, float("-inf"), float("inf")), "Infinite interval")
        self.assertFalse(l2.above_interval(l1, float("-inf"), float("inf")), "Infinite interval")

        #Should intersect on boundary and not be above.
        self.assertFalse(l1.above_interval(l2, float("-inf"), x_intercept),
                         "Intersect on boundary negative infinite test")
        self.assertFalse(l1.above_interval(l2, x_intercept, float("inf")),
                         "Intersect on boundary and below infinite test")
        self.assertFalse(l2.above_interval(l1, float("-inf"), x_intercept),
                         "Intersect on boundary negative infinite test")
        self.assertFalse(l2.above_interval(l1, x_intercept, float("inf")),
                         "Intersect on boundary and below infinite test")

        self.assertTrue(l1.above_interval(l2, float("-inf"), x_intercept - .0001), "Not intersect on boundary and neg inf")
        self.assertFalse(l1.above_interval(l2, x_intercept + .0001, float("inf")), "Not intersect and below infinite")

        self.assertTrue(l2.above_interval(l1, x_intercept + .0001, float("inf")), "Not intersect and above inf")
        self.assertFalse(l2.above_interval(l1, float("-inf"), x_intercept - .0001), "Not intersect and below inf")

# class HLineTestingMethod(unittest.TestCase):
#
#     def test_above_interval(self):
#         l1 = PolyTree.HLine(**{'a': -0.28002249578781524,
#                      'b': 0.04670219515585938})
#         l2 = PolyTree.HLine(**{'a': -0.08642455382757375,
#                      'b': 0.29438253103340406})
#
#         x_intercept = l1.x_intercept(l2)
#
#         self.assertFalse(l1.above_interval(l2, float("-inf"), float("inf")), "Infinite interval")
#         self.assertFalse(l2.above_interval(l1, float("-inf"), float("inf")), "Infinite interval")
#
#         #Should intersect on boundary and not be above.
#         self.assertFalse(l1.above_interval(l2, float("-inf"), x_intercept),
#                          "Intersect on boundary negative infinite test")
#         self.assertFalse(l1.above_interval(l2, x_intercept, float("inf")),
#                          "Intersect on boundary and below infinite test")
#         self.assertFalse(l2.above_interval(l1, float("-inf"), x_intercept),
#                          "Intersect on boundary negative infinite test")
#         self.assertFalse(l2.above_interval(l1, x_intercept, float("inf")),
#                          "Intersect on boundary and below infinite test")
#
#         self.assertTrue(l1.above_interval(l2, float("-inf"), x_intercept - .0001), "Not intersect on boundary and neg inf")
#         self.assertFalse(l1.above_interval(l2, x_intercept + .0001, float("inf")), "Not intersect and below infinite")
#
#         self.assertTrue(l2.above_interval(l1, x_intercept + .0001, float("inf")), "Not intersect and above inf")
#         self.assertFalse(l2.above_interval(l1, float("-inf"), x_intercept - .0001), "Not intersect and below inf")

class SeidelTesting(unittest.TestCase):

        def test_losing_points_cells(self):
            pts = [(random.random(), random.random()) for t in range(1000)]
            tree = compute_cutting(testing_set, {t: 1 for t in testing_set}, pts, 4)
            total_points = set()
            for trap in tree.get_leaves():
                total_points.update(trap.points)
            self.assertEqual(len(total_points), len(pts), "Losing or gaining points somehow")

            t_count = 0
            for trap in tree.get_leaves():
                t_count += len(trap.points)
            self.assertEqual(len(total_points), t_count, "Losing or gaining points somehow")

            all_pts = set(pts)
            for trap in tree.get_leaves():
                all_pts -= set(trap.points)
            self.assertEqual(len(all_pts), 0, "Some of the points have gone missing")

        def test_cell_weights(self):
            """
            Check to see if this is a valid cutting.
            :return:
            """
            pts = [(random.random(), random.random()) for t in range(1000)]
            weight_map = {t: 1 for t in testing_set}
            tree = compute_cutting(testing_set, weight_map, pts, 4)
            total_weight = sum(weight_map[l] for l in testing_set)

            for trap in tree.get_leaves():
                self.assertTrue(trap.get_weight() <= total_weight / 4, "%f <= %f" % (trap.get_weight(), total_weight / 4))

        def test_losing_points_cells_deg(self):
            pts = [(random.random(), random.random()) for t in range(1000)]

            test_set = []
            rand_pts = random.sample(pts, 20)
            for p1, p2 in itertools.combinations(rand_pts, 2):
                test_set.append(to_line(p1, p2))
            random.shuffle(test_set)
            weight_map = {line: 1 for line in test_set}
            tree = compute_cutting(test_set, weight_map, pts, 4)
            total_points = set()
            for trap in tree.get_leaves():
                total_points.update(trap.points)
            self.assertEqual(len(total_points), len(pts), "Losing or gaining points somehow")

            t_count = 0
            for trap in tree.get_leaves():
                t_count += len(trap.points)
            self.assertEqual(len(total_points), t_count, "Losing or gaining points somehow")

            all_pts = set(pts)
            for trap in tree.get_leaves():
                all_pts -= set(trap.points)
            self.assertEqual(len(all_pts), 0, "Some of the points have gone missing")

        def test_cell_weights_deg(self):
            """
            Check to see if this is a valid cutting.
            :return:
            """
            pts = [(random.random(), random.random()) for t in range(1000)]

            test_set = []
            rand_pts = random.sample(pts, 20)
            for p1, p2 in itertools.combinations(rand_pts, 2):
                test_set.append(to_line(p1, p2))
            random.shuffle(test_set)
            weight_map = {line: 1 for line in test_set}
            tree = compute_cutting(test_set, weight_map, pts, 4)
            total_weight = sum(weight_map[l] for l in test_set)


            for trap in tree.get_leaves():
                self.assertTrue(trap.get_weight() <= total_weight / 4,
                                "%f <= %f, \n %s" % (trap.get_weight(), total_weight / 4, str(trap)))

            for trap in tree.get_leaves():
                #print(trap.top_line.x_intercept(trap.bottom_line))
                self.assertTrue(approx_eq_above(trap.bottom_line.evaluate(trap.left_x),
                                                trap.top_line.evaluate(trap.left_x)), "Bad trap left %s"%(str(trap),))
                self.assertTrue(approx_eq_above(trap.bottom_line.evaluate(trap.right_x),
                                                trap.top_line.evaluate(trap.right_x)), "Bad trap right %s"%(str(trap),))

        # def test_cell_incidence(self):
        #     """
        #     Check to see if this is a valid cutting.
        #     :return:
        #     """
        #     pts = [(random.random(), random.random()) for t in range(1000)]
        #
        #     test_set = []
        #     rand_pts = random.sample(pts, 20)
        #     for p1, p2 in itertools.combinations(rand_pts, 2):
        #         test_set.append(to_line(p1, p2))
        #
        #     weight_map = {line: 1 for line in test_set}
        #     tree = compute_cutting(test_set, weight_map, pts, 4)
        #     total_weight = sum(weight_map[l] for l in test_set)
        #
        #
        #     for trap in tree.get_trapezoids():
        #         self.assertTrue(trap.getWeight() <= total_weight / 4,
        #                         "%f <= %f, \n %s" % (trap.getWeight(), total_weight / 4, str(trap)))
        #
        #     for trap in tree.get_trapezoids():
        #         #print(trap.top_line.x_intercept(trap.bottom_line))
        #         self.assertTrue(approx_eq_above(trap.bottom_line.evaluate(trap.left_x),
        #                                         trap.top_line.evaluate(trap.left_x)), "Bad trap left %s"%(str(trap),))
        #         self.assertTrue(approx_eq_above(trap.bottom_line.evaluate(trap.right_x),
        #                                         trap.top_line.evaluate(trap.right_x)), "Bad trap right %s"%(str(trap),))


class ConeTesting(unittest.TestCase):
    pass

class PolyTesting(unittest.TestCase):

    def setUp(self):
        self.pts = [(random.random(), random.random()) for t in range(1000)]

        self.test_set = []
        rand_pts = random.sample(self.pts, 40)
        for p1, p2 in itertools.combinations(rand_pts, 2):
            self.test_set.append(to_line(p1, p2))
        random.shuffle(self.test_set)
        self.weight_map = {line: 1 for line in self.test_set}
        self.r = 9

        self.tree = poly_tree.compute_cutting(self.test_set, self.weight_map, self.pts, self.r)
        self.total_weight = sum(self.weight_map[l] for l in self.test_set)

    def test_losing_points_cells(self):
        total_points = set()
        for trap in self.tree.get_leaves():
            total_points.update(trap.points)
        self.assertEqual(len(total_points), len(self.pts), "Losing or gaining points somehow")

        t_count = 0
        for trap in self.tree.get_leaves():
            t_count += len(trap.points)
        self.assertEqual(len(total_points), t_count, "Losing or gaining points somehow")

        all_pts = set(self.pts)
        for trap in self.tree.get_leaves():
            all_pts -= set(trap.points)
        self.assertEqual(len(all_pts), 0, "Some of the points have gone missing")

    def test_cell_weights(self):
        """
        Check to see if this is a valid cutting.
        :return:
        """
        for trap in self.tree.get_leaves():
            self.assertFalse(trap.total_weight(self.weight_map) > sum(self.weight_map[l] for l in self.weight_map) / self.r,
                            "%f <= %f" % (trap.total_weight(self.weight_map), self.total_weight / self.r))

    def test_losing_points_cells_deg(self):

        total_points = set()
        for trap in self.tree.get_leaves():
            total_points.update(trap.points)
        self.assertEqual(len(total_points), len(self.pts), "Losing or gaining points somehow")

        t_count = 0
        for trap in self.tree.get_leaves():
            t_count += len(trap.points)
        self.assertEqual(len(total_points), t_count, "Losing or gaining points somehow")

        all_pts = set(self.pts)
        for trap in self.tree.get_leaves():
            all_pts -= set(trap.points)
        self.assertEqual(len(all_pts), 0, "Some of the points have gone missing")

    def test_wedge(self):
        dual_set = [to_dual_line(pt) for pt in self.pts]
        x1 = random.random()
        x2 = random.random()


        query_line = Segment(Line(.2 * random.random() + .8 , -random.random()), min(x1, x2), max(x1, x2))
        query_wedge = l_wedge(query_line)

        overlap_set = []
        for c in self.tree.get_leaves():
            for p in c.get_border_vertices():
                if math.isfinite(p[0]) and query_wedge.contains_pt(p):
                    overlap_set.append(p)
        overlap_set = deduplicate_points(overlap_set)
        overlap_set_tree = self.tree.count_wedge(query_wedge)
        overlap_set_tree = deduplicate_points(overlap_set_tree)


        f, ax = plt.subplots()
        self.tree.visualize_arrangement(ax, -1, 2, -1, 2)

        if len(overlap_set) > 0:
            x, y = zip(*overlap_set)
            ax.scatter(x, y, color='b')
        for c in self.tree.wedge_projection(query_wedge):
            c.visualize(ax, -1, 2, -1, 2, color='r')

        lx = min(max(query_wedge.up_segment.xl, -1), 2)
        rx = min(max(query_wedge.up_segment.xr, -1), 2)

        ax.plot([lx, rx], [query_wedge.up_segment.evaluate(lx), query_wedge.up_segment.evaluate(rx)], c='r')
        ax.plot([lx, rx], [query_wedge.down_segment.evaluate(lx), query_wedge.down_segment.evaluate(rx)], c='b')
        ax.set_xlim([-1, 2])
        ax.set_ylim([-1, 2])
        plt.show()

        self.assertEqual(len(overlap_set), len(overlap_set_tree))








# def test_trap_case():
#
#     matplotlib.rcParams['figure.figsize'] = [30.0, 30.0]
#     r = 4
#     pts = [(random.random(), random.random()) for t in range(1000)]
#
#     w_lines = []
#     rand_pts = random.sample(pts, 20)
#     for p1, p2 in itertools.combinations(rand_pts, 2):
#         w_lines.append((to_line(p1, p2), 1))
#
#
#     min_weight = len(w_lines) / float(r)
#     tree = Seidel_Tree(w_lines, points=pts, min_weight=min_weight, min_p_count=-1)
#     for l, _ in w_lines:
#         f, a = plt.subplots()
#         tree.insert_line(l)
#         tree.visualize_trapezoids(a, 0, .9, -1, 1)
#         plt.show()
#
# test_trap_case()

# def viz_degenerate_case():
#
#     pts = [(random.random(), random.random()) for t in range(10)]
#     testing_set = []
#     for p in pts[1:4]:
#         testing_set.append(to_line(pts[0], p))
#     for p in pts[6:]:
#         testing_set.append(to_line(p, pts[5]))
#     tree = Seidel_Tree([], [], min_weight=-1)
#     for l in testing_set:
#         tree.insert_line(l)
#     tree.visualize("output_set")
#     matplotlib.rcParams['figure.figsize'] = [30.0, 30.0]
#     f, a = plt.subplots()
#     tree.visualize_arrangement(a, -1, 1, -1, 1)
#     plt.show()
#
# viz_degenerate_case()


# class CellTestingMethod(unittest.TestCase):
#     """
#     Need to test to make sure all the lines added overlap the correct cells.
#     """
#
#     def test_insert_splitter(self):


if __name__ == '__main__':
    unittest.main()
