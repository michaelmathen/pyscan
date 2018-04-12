import random
import SeidelTree
import PolyTree
import test_set
import line_testing
import Partitioning
import csv
import numpy as np
import time
import math


# Plots showing sample size vs. error
# Plots showing sample size vs. time
# Plots showing test set size vs. error
# Plots showing r vs. error
# Plots showing test set size vs. time
# Plots showing r vs. time


def testing_framework(pts, l_s, h_s, count, output_file = None,
                      vparam="output_size", r=4, input_size=None,
                      output_size=200, cell_size=1,
                      test_set_f="lts",
                      cutting_f="poly",
                      part_f="pchan",
                      k=100):

    if output_file is None:
        if part_f == "box" or part_f == "sample":
            output_file = vparam + "_" + part_f
        else:
            output_file = vparam + "_" + part_f + "_" + cutting_f + "_" + test_set_f + "_"
            if vparam != "r":
                output_file += "_" + str(r)
        if vparam != input_size:
            if input_size is None:
                output_file += "_" + str(len(pts))
            else:
                output_file += "_" + str(input_size)
        if vparam != output_size:
            output_file += "_" + str(output_size)
    if vparam != cell_size and part_f != "sample":
        output_file += "_" + str(cell_size)

    fieldnames = ["vparam", "r", "input_size", "output_size", "cell_size",
                  "test_set_f", "cutting_f", "part_f", "time", "error", "k"]
    if input_size is None:
        input_size = len(pts)
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i in np.linspace(l_s, h_s, count):
            if vparam == "r":
                r = max(int(i + .5), 2)
            elif vparam == "input_size":
                input_size = max(int(i + .5), 1)
            elif vparam == "output_size":
                output_size = max(int(i + .5), 1)
            elif vparam == "cell_size":
                cell_size =  max(int(i + .5), 1)
            else:
                raise ValueError("vparam=%s"%(vparam,))

            if cutting_f == "poly":
                compute_cutting = PolyTree.compute_cutting
            elif cutting_f == "trap":
                compute_cutting = SeidelTree.compute_cutting
            else:
                raise ValueError("cutting_f=%s"%(cutting_f, ))

            if test_set_f == "lts":
                test_f = test_set.test_set_lines
            elif test_set_f == "pts":
                test_f = test_set.test_set_points
            elif test_set_f == "dts":
                test_f = test_set.test_set_dual
            elif test_set_f == "dts_x":
                test_f = test_set.test_set_dual_exact_t
            else:
                raise ValueError("test_set_f=%s"%(test_set_f,))

            input_set = random.sample(pts, input_size)
            min_cell_size = input_size / output_size * cell_size
            #(min_cell_size, len(input_set), r, cell_size)
            start_time = time.time()
            if part_f == "pchan":
                output_set, weights = Partitioning.chan_partitions(input_set, r,
                                             min_cell_size=min_cell_size,
                                             cell_sample_size=cell_size,
                                             cutting_f=compute_cutting,
                                             test_set_f=test_f, c=.05)
            elif part_f == "pmat":
                output_set, weights = Partitioning.partitions(input_set, r,
                                                     min_cell_size=min_cell_size,
                                                     cell_sample_size=cell_size,
                                                     cutting_f=compute_cutting,
                                                     test_set_f=test_f)
            elif part_f == "sample":
                output_set = random.sample(input_set, output_size)
                weights = [1.0 for p in output_set]
            elif part_f == "box":
                output_set, weights = Partitioning.quadTreeSample(input_set, min_cell_size)
            else:
                raise ValueError("part_f=%s"%(part_f,))
            end_time = time.time()

            max_error = line_testing.line_error_test(pts, output_set, weights)

            row = {"vparam": vparam,
                   "r" : r,
                             "input_size": input_size,
                             "output_size":len(output_set),
                             "cell_size": min_cell_size,
                             "test_set_f": test_set_f,
                             "cutting_f": cutting_f,
                             "part_f": part_f,
                             "time": end_time - start_time,
                             "error": max_error,
                             "k": k}
            writer.writerow(row)
            print(row)



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

if __name__ == "__main__":
    pts = [(random.random(), random.random()) for i in range(100000)]
    #pts = upload_crimes("crimes.csv")
    print(len(pts))
    #Vary the output sample size
    for cutting in ["poly"]:
        for l_name in [ "dts_x"]:
            testing_framework(pts, 50, 1000, 10, part_f="pchan",
                              test_set_f=l_name,
                              cutting_f=cutting,
                              r=2)
    #
    # #Vary the r parameter used
    # for cutting in ["poly", "trap"]:
    #     for l_name in ["lts", "pts"]:
    #         testing_framework(pts, 2, 12, 10, vparam="r", part_f="pchan",
    #                           test_set_f=l_name, cutting_f=cutting)

    #Vary the input size with fixed output size
    # for cutting in ["poly"]:
    #     for l_name in ["lts", "pts"]:
    #         testing_framework(pts, 1000, 100000, 10, vparam="input_size", part_f="pchan",
    #                           test_set_f=l_name, cutting_f=cutting,
    #                           r=2)



    # testing_framework(pts, 50, 2000, 10, part_f="box")
    # testing_framework(pts, 1000, 100000, 10, vparam="input_size", part_f="box")
    #
    # #
    # testing_framework(pts, 200, 2000, 10, part_f="sample", output_file="output_sampling_chan.csv")
    # #
    # testing_framework(pts, 1000, 100000, 10, part_f="sample", vparam="input_size", output_file="input_sampling_chan.csv")


    print("Done")
