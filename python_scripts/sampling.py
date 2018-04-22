import random
import numpy as np
import math

"""
File containing sampling related codes.
"""


def weighted_shuffle(items, weights):
    rweights = [-random.random() ** (1.0 / w) for w in weights]
    order = sorted(range(len(items)), key=lambda i: rweights[i])
    return [items[i] for i in order]


def random_gap(p):
    u = max(random.random(), np.finfo(float).eps)
    return math.floor(math.log(u) / math.log(1-p))


def p_sample(data, p):
    if p > .3:
        for el in data:
            if random.random() < p:
                yield el
    else:
        i = random_gap(p)
        while i < len(data):
            yield data[i]
            i += random_gap(p)


def binom_rand_k(n, p):
    return sum(p_sample([1] * int(round(n) +.1), p))
