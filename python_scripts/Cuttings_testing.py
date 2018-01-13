import Cuttings
import unittest
import matplotlib.pyplot as plt



class LineTestingMethod(unittest.TestCase):

    def test_above_interval(self):
        l1 = Cuttings.Line(**{'a': -0.28002249578781524,
                     'b': 0.04670219515585938})
        l2 = Cuttings.Line(**{'a': -0.08642455382757375,
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

# class CellTestingMethod(unittest.TestCase):
#     """
#     Need to test to make sure all the lines added overlap the correct cells.
#     """
#
#     def test_insert_splitter(self):


if __name__ == '__main__':
    unittest.main()
