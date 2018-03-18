import random
from SeidelTree import to_line


def line_error_test(big_samp, small_samp, weights, pt_num):
    test_pool = big_samp + small_samp
    for i in range(pt_num):
        [p1, p2] = random.sample(test_pool, 2)
        try:
            line = to_line(p1, p2)
            r_1 = sum(line.pt_eq_below(p) for p in big_samp) / float(len(big_samp))
            r_2 = sum(weights[i] for i, p in enumerate(small_samp) if line.pt_eq_below(p)) / sum(weights[i] for i, p in enumerate(small_samp))
        except ZeroDivisionError:
            continue

        yield abs(r_1 - r_2)


def line_error_test_unweighted(big_samp, small_samp, pt_num):
    test_pool = big_samp + small_samp
    for i in range(pt_num):
        [p1, p2] = random.sample(test_pool, 2)
        try:
            line = to_line(p1, p2)
            r_1 = sum(line.pt_eq_below(p) for p in big_samp) / float(len(big_samp))
            r_2 = sum(line.pt_eq_below(p) for p in small_samp) / float(len(small_samp))
        except ZeroDivisionError:
            continue

        yield abs(r_1 - r_2)
