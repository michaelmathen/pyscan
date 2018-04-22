
import test_set



def biased_l2(pts, alpha):
    r = int(round(alpha * len(pts)) + .1)
    lines = test_set.test_set_points(pts, r * r)
    sample = []
    range_sample_counts = {l:0  for l in lines}
    range_data_counts = {l:0 for l in lines}

    pt_sample_counts = {pt:0 for pt in pts}
    pt_data_counts = {pt:0 for pt in pts}
    for p in pts:
        size = 0
        for l in lines:
            if l.pt_eq_below(p):
                range_data_counts[l] += 1
                pt_sample_counts[p] += range_sample_counts[l]
                pt_data_counts[p] += range_data_counts[l]
                size += 1

        if size / 2 + pt_sample_counts[p] - alpha * pt_data_counts[p] <= 0:
            sample.append(p)
            for l in lines:
                if l.pt_eq_below(p):
                    range_sample_counts[l] += 1

    return sample


