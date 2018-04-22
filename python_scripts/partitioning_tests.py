import random
import seidel_tree
import poly_tree
import test_set
import line_testing
import partitioning
import csv
import numpy as np
import time
import math
import akcan


# Plots showing sample size vs. error
# Plots showing sample size vs. time
# Plots showing test set size vs. error
# Plots showing r vs. error
# Plots showing test set size vs. time
# Plots showing r vs. time


def testing_framework(pts, l_s, h_s, count, output_file = None,
                      vparam="output_size", b=32, input_size=None,
                      output_size=500, cell_size=1,
                      test_set_f="lts",
                      cutting_f="poly",
                      part_f="pchan"):

    if output_file is None:
        if part_f == "box" or part_f == "sample" or part_f == "l2":
            output_file = vparam + "_" + part_f + ".csv"
        else:
            output_file = vparam + "_" + part_f + "_" + cutting_f + "_" + test_set_f + ".csv"


    fieldnames = ["vparam", "b", "input_size", "output_size", "cell_size",
                  "test_set_f", "cutting_f", "part_f", "time", "error", "k"]
    if input_size is None:
        input_size = len(pts)
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i in np.linspace(l_s, h_s, count):
            if vparam == "b":
                b = max(int(i + .5), 2)
            elif vparam == "input_size":
                input_size = max(int(i + .5), 1)
            elif vparam == "output_size":
                output_size = max(int(i + .5), 1)
            elif vparam == "cell_size":
                cell_size =  max(int(i + .5), 1)
            else:
                raise ValueError("vparam=%s"%(vparam,))

            input_set = random.sample(pts, input_size)
            min_cell_size = input_size / output_size * cell_size
            #(min_cell_size, len(input_set), r, cell_size)
            start_time = time.time()
            if part_f == "pchan":
                tree = partitioning.chan_partitions2(input_set, b,
                                                     min_cell_size=min_cell_size)
                output_set, weights = partitioning.sample_cells([leaf.get_points() for leaf in tree.get_leaves()],
                                                                cell_size)
            elif part_f == "pchans":
                tree = partitioning.chan_partitions_simple(input_set, b,
                                                     min_cell_size=min_cell_size)
                output_set, weights = partitioning.sample_cells([leaf.get_points() for leaf in tree.get_leaves()],
                                                                cell_size)
            elif part_f == "pmat":
                # if cutting_f == "poly":
                #     compute_cutting = poly_tree.compute_cutting
                # elif cutting_f == "trap":
                #     compute_cutting = seidel_tree.compute_cutting
                # else:
                #     raise ValueError("cutting_f=%s" % (cutting_f,))

                if test_set_f == "lts":
                    test_f = test_set.test_set_lines
                elif test_set_f == "pts":
                    test_f = test_set.test_set_points
                elif test_set_f == "dts":
                    test_f = test_set.test_set_dual_exact_t
                else:
                    raise ValueError("test_set_f=%s" % (test_set_f,))
                output_set, weights = partitioning.partitions(input_set, b,
                                                              min_cell_size=min_cell_size,
                                                              cell_sample_size=cell_size,
                                                              test_set_f=test_f)
            elif part_f == "sample":
                output_set = random.sample(input_set, output_size)
                weights = [1.0 for _ in output_set]
            elif part_f == "box":
                output_set, weights = partitioning.quadTreeSample(input_set, min_cell_size)
            elif part_f == "l2":
                output_set = akcan.biased_l2(input_set, output_size / input_size)
                weights = [1.0 for _ in output_set]
            else:
                raise ValueError("part_f=%s"%(part_f,))
            end_time = time.time()

            max_error = line_testing.line_error_test(pts, output_set, weights)
            row = {"vparam": vparam,
                   "b" : b,
                     "input_size": input_size,
                     "output_size":len(output_set),
                     "cell_size": min_cell_size,
                     "test_set_f": test_set_f,
                     "cutting_f": cutting_f,
                     "part_f": part_f,
                     "time": end_time - start_time,
                     "error": max_error}
            writer.writerow(row)
            print(row)
    return len(output_set)

def upload_crimes(file_name, perturb = 1):
    with open(file_name) as f:
        reader = csv.DictReader(f)
        pts = []
        theta = random.random() * math.pi

        a1 = math.cos(theta)
        a2 = -math.sin(theta)
        b1 = math.sin(theta)
        b2 = math.cos(theta)
        for row in reader:
            try:
                x = float(row["X Coordinate"]) + random.random() * perturb
                y = float(row["Y Coordinate"]) + random.random() * perturb
                pts.append((a1 * x + a2 * y, b1 * x + b2 * y))
            except ValueError:
                pass
        return pts

def b_test(pts, b_low, b_high, output_size, count=20):
    sample_size = testing_framework(pts, b_low, b_high, count, part_f="pchans", vparam="b",
                      test_set_f="dts",
                      cutting_f="poly", output_size=output_size)
    sample_size = testing_framework(pts, b_low, b_high, count, part_f="pchan", vparam="b",
                      test_set_f="dts",
                      cutting_f="poly", output_size=output_size)
    #
    # for cutting in ["poly"]:
    #     for l_name in ["dts", "lts", "pts"]:
    #         testing_framework(pts, b_low, b_high, count, part_f="pmat", vparam="b",
    #                           test_set_f=l_name,
    #                           cutting_f=cutting, output_size=output_size)


def output_size(pts, output_low, output_high, count=10):
    #testing output size
    testing_framework(pts, output_low, output_high, count, b=22, part_f="pchan",
                        test_set_f="dts",
                        cutting_f="poly")
    testing_framework(pts, output_low, output_high, count, b=22, part_f="pchans",
                                    test_set_f="dts",
                                    cutting_f="poly")
    #testing_framework(pts, output_low, output_high, count, part_f="sample")
    # testing_framework(pts, output_low, output_high, count, part_f="box")
    #
    # for l_name in [ "dts", "lts", "pts"]:
    #     testing_framework(pts, output_low, output_high, count, part_f="pmat",
    #                       test_set_f=l_name,
    #                       cutting_f="poly",
    #                       b=16)

def input_size(pts, input_low, input_high, output_size, count=10):
    #Vary the input size with fixed output size
    # for cutting in ["poly"]:
    #     for l_name in ["dts", "lts", "pts"]:
    #         testing_framework(pts, input_low, input_high, count, vparam="input_size", part_f="pmat",
    #                           test_set_f=l_name, cutting_f=cutting, output_size=output_size,
    #                           b=16)

    sample_size = testing_framework(pts, input_low, input_high, count, b=22, part_f="pchan", output_size=output_size, vparam="input_size")
    sample_size = testing_framework(pts, input_low, input_high, count, b=22, part_f="pchans", output_size=output_size, vparam="input_size")
    #
    # testing_framework(pts, input_low, input_high, count, part_f="sample", output_size=output_size, vparam="input_size")
    # testing_framework(pts, input_low, input_high, count, part_f="box", vparam="input_size")


if __name__ == "__main__":
    #pts = [(random.random(), random.random()) for i in range(100000)]
    pts = upload_crimes("crimes.csv")
    default_output = 500

    pts = random.sample(pts, 100000)
    #
    #b_test(pts, 4, 32, default_output)
    output_size(pts, 25, 2000)
    input_size(pts, 10000, len(pts), output_size=default_output)

    print("Done")
