import random
import math
import bisect
from SeidelTree import to_line


def line_error_test(big_samp, small_samp, weights, n=None):
    return line_error_test_fast(big_samp, small_samp, weights, n)


def line_error_test_unweighted(big_samp, small_samp, n=None):
    return line_error_test_fast(big_samp, small_samp, [1 for p in small_samp], n)


def order_function(p1, p2):
    x = p2[0] - p1[0]
    y = p2[1] - p1[1]
    if y >= 0:
        return math.atan2(y, x)
    else:
        return math.pi - math.atan2(y, x)


def line_error_test_fast(big_set, small_set, weights, n=None):
    """
    :param big_samp:
    :param small_samp:
    :param weights:
    :param pt_num:
    :param n: The sub-sampling size of the small sample.
    :return:
    """

    if n is not None:
        net_sample = random.sample(small_set, n)
    else:
        net_sample = small_set[:]

    for i in range(len(net_sample)):
        sample_part = net_sample[i + 1:]

        order_f = lambda x: order_function(net_sample[i], x)

        deltas_big = [0] * (len(sample_part) + 1)
        deltas_small = [0] * (len(sample_part) + 1)
        angles = [order_f(p) for p in sample_part]
        angles.sort()

        for p_1 in big_set:
            insertion_pt = bisect.bisect_left(angles, order_f(p_1))
            deltas_big[insertion_pt] += math.copysign(1, p_1[1])

        for p_1, w in zip(small_set, weights):
            insertion_pt = bisect.bisect_left(angles, order_f(p_1))
            deltas_small[insertion_pt] += math.copysign(w, p_1[1])

        big_sample_curr = sum(1 for p in big_set if p[1] <= 0)
        small_sample_curr = sum(w for p, w in zip(small_set, weights) if p[1] <= 0)

        max_discrepancy = 0
        a = 1.0 / len(big_set)
        b = 1.0 / sum(weights)
        for db, ds in zip(deltas_big, deltas_small):
            big_sample_curr += db
            small_sample_curr += ds
            max_discrepancy = max(abs(big_sample_curr * a - small_sample_curr * b), max_discrepancy)
        return max_discrepancy




