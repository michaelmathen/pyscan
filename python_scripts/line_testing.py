import random
import math
import bisect
from seidel_tree import to_line
import line_discrepancy

def line_error_test(big_samp, small_samp, weights, n=None):
    return line_error_test_fast(big_samp, small_samp, weights, n)


def line_error_test_unweighted(big_samp, small_samp, n=None):
    return line_error_test_fast(big_samp, small_samp, [1 for p in small_samp], n)


def line_error_test_fast(big_set, small_set, weights, n=None):
    """
    :param big_samp:
    :param small_samp:
    :param weights:
    :param pt_num:
    :param n: The sub-sampling size of the small sample.
    :return:
    """
    net = random.sample(big_set, len(small_set) // 4)
    md, ml = line_discrepancy.line_discrepancy(net, big_set, [1 for _ in big_set], small_set, weights, lambda m, b: abs(m - b))
    return md



