import random
import SeidelTree
import PolyTree
import test_set
import line_testing
import Partitioning
import csv
import numpy as np
import time


# Plots showing sample size vs. error
# Plots showing sample size vs. time
# Plots showing test set size vs. error
# Plots showing r vs. error
# Plots showing test set size vs. time
# Plots showing r vs. time


def testing_framework(pts, l_s, h_s, count, output_file,
                      vparam="output_size", r=4, input_size=10000,
                      output_size=200, cell_size=1,
                      test_set_f="lts",
                      cutting_f="poly",
                      part_f="pchan",
                      k=100):

    fieldnames = ["vparam", "r", "input_size", "output_size", "cell_size",
                  "test_set_f", "cutting_f", "part_f", "time", "error", "k"]
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
                                             test_set_f=test_f)
            elif part_f == "pmat":
                output_set, weights = Partitioning.partitions(input_set, r,
                                                     c = .1,
                                                     min_cell_size=min_cell_size,
                                                     cell_sample_size=cell_size,
                                                     cutting_f=compute_cutting,
                                                     test_set_f=test_f)
            elif part_f == "sample":
                output_set = random.sample(input_set, output_size)
                weights = [1.0 for p in output_set]
            else:
                raise ValueError("part_f=%s"%(part_f,))
            end_time = time.time()

            max_error = 0
            for e in line_testing.line_error_test(pts, output_set, weights, k):
                max_error = max(max_error, e)
            row = {"vparam": vparam,
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


if __name__ == "__main__":
    pts = [(random.random(), random.random()) for i in range(100000)]
    # for cutting in ["poly", "trap"]:
    #     for l_name in ["lts", "pts"]:
    #         testing_framework(pts, 50, 1000, 10, "timing_plot_r2_"+ cutting + "_" + l_name + ".csv",
    #                           test_set_f=l_name, cutting_f=cutting, r=2)
    #
    # for cutting in ["poly", "trap"]:
    #     for l_name in ["lts", "pts"]:
    #         testing_framework(pts, 50, 1000, 10, "timing_plot_r4_"+ cutting + "_" + l_name + ".csv",
    #                           test_set_f=l_name, cutting_f=cutting, r=4)

    testing_framework(pts, 200, 1800, 10, "sampling.csv", part_f="sample")

    print("Done")
#
# def random_sample_plot(point_count, l_s, h_s, k):
#     big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
#     errors = []
#     sizes = []
#     for i in np.linspace(l_s, h_s, k):
#         s_size = int(i + .5)
#         sampled_pts = random.sample(big_sample, s_size)
#         #print(sampled_weights)
#         max_error = 0
#         for e in line_error_test_unweighted(big_sample, sampled_pts,  100):
#             max_error = max(max_error, e)
#         print("finished run")
#         errors.append(max_error)
#         sizes.append(len(sampled_pts))
#     f, a = plt.subplots()
#     a.scatter(sizes, errors)
#     plt.show()
#
# def random_sample_time_plot(point_count, l_s, h_s, k):
#     big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
#     times = []
#     sizes = []
#     for i in np.linspace(l_s, h_s, k):
#         s_size = int(i + .5)
#         s = time()
#         sampled_pts = random.sample(big_sample, s_size)
#         e = time()
#         times.append(e - s)
#         sizes.append(len(sampled_pts))
#     f, a = plt.subplots()
#     a.scatter(sizes, times)
#     plt.show()
#
# def random_sample_time_plot2(point_count, l_s, h_s, k):
#     times = []
#     sizes = []
#     for i in np.linspace(l_s, h_s, k):
#         big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(int(i))]
#         s = time()
#         sampled_pts = random.sample(big_sample, point_count)
#         e = time()
#         times.append(e - s)
#         sizes.append(i)
#     f, a = plt.subplots()
#     a.scatter(sizes, times)
#     plt.show()
#
#
# def error_plot(point_count, l_s, h_s, k, r=4, c=2):
#     big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
#     errors = []
#     sizes = []
#     for i in np.linspace(l_s, h_s, k):
#         sampled_pts, sampled_weights = partitions(big_sample, r, c, min_cell_size=i)
#         #print(sampled_weights)
#         max_error = 0
#         for e in line_error_test(big_sample, sampled_pts, sampled_weights, 100):
#             max_error = max(max_error, e)
#         print("finished run %f"%(i, ))
#         errors.append(max_error)
#         sizes.append(len(sampled_pts))
#     #print(sizes, errors)
#     f, a = plt.subplots()
#     a.scatter(sizes, errors)
#     plt.show()
#
#
# def time_plot(point_count, min_cell_size, l_s, h_s, k, r=4, c=2):
#     big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
#     times = []
#     sizes = []
#     for i in np.linspace(l_s, h_s, k):
#         s = time()
#         sampled_pts, sampled_weights = partitions(big_sample, r, c, min_cell_size=min_cell_size, cell_sample_size=i,
#                                                   test_set_f=random_test_set)
#         e = time()
#         #print(sampled_weights)
#         print("finished run %f"%(i, ))
#         times.append(e - s)
#         sizes.append(len(sampled_pts))
#     #print(sizes, errors)
#     f, a = plt.subplots()
#     a.scatter(sizes, times)
#     plt.show()
#
# def time_plot2(output_size, l_s, h_s, k, r=4, c=2):
#     times = []
#     sizes = []
#     for i in np.linspace(l_s, h_s, k):
#         big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(int(i))]
#         min_cell_size = i / output_size
#         s = time()
#         sampled_pts, sampled_weights = chan_partitions(big_sample, r,
#                                                   min_cell_size=min_cell_size,
#                                                   cell_sample_size=1,
#                                                     cutting_f=Pt.compute_cutting,
#                                                   test_set_f=random_test_set)
#         e = time()
#         #print(sampled_weights)
#         print("finished run %f"%(i, ))
#         times.append(e - s)
#         sizes.append(i)
#     #print(sizes, errors)
#     f, a = plt.subplots()
#     a.scatter(sizes, times)
#     plt.show()
#
# def error_plot(output_size, l_s, h_s, k, r=4):
#     sizes = []
#     errors = []
#     for i in np.linspace(l_s, h_s, k):
#         big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(int(i))]
#         min_cell_size = i / output_size
#         sampled_pts, sampled_weights = chan_partitions(big_sample, r,
#                                                   min_cell_size=min_cell_size,
#                                                   cell_sample_size=1,
#                                                     cutting_f=Pt.compute_cutting,
#                                                   test_set_f=random_test_set)
#         max_error = 0
#         for e in line_error_test(big_sample, sampled_pts, sampled_weights, 100):
#             max_error = max(max_error, e)
#         errors.append(max_error)
#         sizes.append(i)
#     #print(sizes, errors)
#     f, a = plt.subplots()
#     a.scatter(sizes, errors)
#     plt.show()
#
#
# def error_plot_output_size(set_size, l_s, h_s, k, r=4):
#     sizes = []
#     errors = []
#     big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(set_size)]
#     for i in np.linspace(l_s, h_s, k):
#
#         min_cell_size = set_size / i
#         sampled_pts, sampled_weights = chan_partitions(big_sample, r,
#                                                   min_cell_size=min_cell_size,
#                                                   cell_sample_size=1,
#                                                     cutting_f=Pt.compute_cutting,
#                                                   test_set_f=random_test_set)
#         max_error = 0
#         for e in line_error_test(big_sample, sampled_pts, sampled_weights, 100):
#             max_error = max(max_error, e)
#         errors.append(max_error)
#         sizes.append(len(sampled_pts))
#     #print(sizes, errors)
#     f, a = plt.subplots()
#     a.scatter(sizes, errors)
#     plt.show()