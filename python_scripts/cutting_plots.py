import random
import numpy as np
import matplotlib.pyplot as plt
import test_set
import PolyTree
import SeidelTree
import PolyTree as pt
import SeidelTree as ts
import time
import test_set
import csv

"""
#377eb8
#4daf4a
#984ea3
"""
cutting_map = {"tsd_color" : "#e41a1c",
               "psd_color" : "b",
               "lsd_color" : "k",
               "poly_marker" : "x",
               "trap_marker" : "o",
               "poly_name": "poly cutting",
                "trap_name": "trap_cutting"}


# def testing_framework(pts, l_s, h_s, count, output_file,
#                       r=2,
#                       line_count=1000,
#                       vparam="line_count",
#                       test_set_f="lts",
#                       cutting_f="poly"):
#
#     fieldnames = ["vparam", "r", "line_count", "output_size",
#                   "test_set_f", "cutting_f", "time"]
#     with open(output_file, 'w') as f:
#         writer = csv.DictWriter(f, fieldnames=fieldnames)
#         writer.writeheader()
#
#         for i in np.linspace(l_s, h_s, count):
#             if vparam == "r":
#                 r = max(int(i + .5), 2)
#             elif vparam == "line_count":
#                 line_count = max(int(i + .5), 1)
#             else:
#                 raise ValueError("vparam=%s"%(vparam,))
#
#
#
#             if test_set_f == "lts":
#                 test_f = test_set.test_set_lines
#             elif test_set_f == "pts":
#                 test_f = test_set.test_set_points
#             elif test_set_f == "dts":
#                 test_f = test_set.test_set_dual
#             else:
#                 raise ValueError("test_set_f=%s"%(test_set_f,))
#
#             start_time_lines = time.time()
#             lines = test_f(pts, line_count)
#             end_time_lines = time.time()
#
#
#             row = {"vparam": vparam,
#                             "line_count": line_count,
#                              "output_size":len(output_set),
#                              "cell_size": min_cell_size,
#                              "test_set_f": test_set_f,
#                              "cutting_f": cutting_f,
#                              "part_f": part_f,
#                              "time": end_time - start_time,
#                              "error": max_error,
#                              "k": k}
#             writer.writerow(row)
#             print(row)

def plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f, test_set_f, r=4, marker='.', color='r'):

    line_count = []
    trap_count = []
    for i in np.linspace(l_s, h_s, k):
        print(i)
        test_set = test_set_f(pts, int(i))
        weight_map = {t: 1 for t in test_set}
        tree = cutting_f(test_set, weight_map, pts, r)
        trap_count.append(len(tree.get_leaves()) / float(r * r))
        line_count.append(len(test_set))

    ax.scatter(line_count, trap_count, marker=marker, c=color)


def plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f, test_set_f, r=4, marker='.', color='r'):

    line_count = []
    times = []
    for i in np.linspace(l_s, h_s, k):
        print(i)
        test_set = test_set_f(pts, int(i))
        weight_map = {t: 1 for t in test_set}
        t = time.time()
        tree = cutting_f(test_set, weight_map, pts, r)
        e = time.time()
        times.append(e - t)
        line_count.append(len(test_set))

    ax.scatter(line_count, times, marker=marker, c=color)


def generate_cutting_size_test(pts, l_s, h_s, k, r=4):
    f, ax = plt.subplots()

    # print("dual")
    # plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=ts.test_set_dual, r=r,
    #                   marker=cutting_map["poly_marker"],
    #                   color=cutting_map["tsd_color"])
    print("points ")
    plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting_greedy, test_set_f=test_set.test_set_points, r=r,
                      marker=cutting_map["poly_marker"],
                      color=cutting_map["psd_color"])
    # print("lines")
    # plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_lines, r=r,
    #                   marker=cutting_map["poly_marker"],
    #                   color=cutting_map["lsd_color"])

    # print("dual")
    # plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=ts.test_set_dual, r=r,
    #                   marker=cutting_map["trap_marker"],
    #                   color=cutting_map["tsd_color"])
    # print("points")
    # plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=ts.test_set_points, r=r,
    #                   marker=cutting_map["trap_marker"],
    #                   color=cutting_map["psd_color"])
    # print("lines")
    # plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=ts.test_set_lines, r=r,
    #                   marker=cutting_map["trap_marker"],
    #                   color=cutting_map["lsd_color"])
    #    ax.legend(["poly dts", "poly psd", "poly lsd", "trap dts", "trap psd", "trap lsd"])
    ax.legend(["poly psd", "poly lsd",  "trap psd", "trap lsd"])
    ax.set_xlabel("Line number")
    ax.set_ylabel("Cutting constant")
    plt.show()


def generate_cutting_time_test(pts, l_s, h_s, k, r=4):
    f, ax = plt.subplots()

    # print("dual")
    # plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=ts.test_set_dual, r=r,
    #                   marker=cutting_map["poly_marker"],
    #                   color=cutting_map["tsd_color"])
    print("points ")
    plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting_greedy, test_set_f=test_set.test_set_points, r=r,
                      marker=cutting_map["trap_marker"],
                      color=cutting_map["tsd_color"])
    print("points ")
    plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_points, r=r,
                      marker=cutting_map["poly_marker"],
                      color=cutting_map["psd_color"])
    # print("lines")
    # plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=ts.test_set_lines, r=r,
    #                   marker=cutting_map["poly_marker"],
    #                   color=cutting_map["lsd_color"])

    # print("dual")
    # plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=ts.test_set_dual, r=r,
    #                   marker=cutting_map["trap_marker"],
    #                   color=cutting_map["tsd_color"])
    # print("points")
    # plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=ts.test_set_points, r=r,
    #                   marker=cutting_map["trap_marker"],
    #                   color=cutting_map["psd_color"])
    # print("lines")
    # plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=ts.test_set_lines, r=r,
    #                   marker=cutting_map["trap_marker"],
    #                   color=cutting_map["lsd_color"])
    #    ax.legend(["poly dts", "poly psd", "poly lsd", "trap dts", "trap psd", "trap lsd"])
    ax.legend(["poly greedy", "poly greedy"])#,  "trap psd", "trap lsd"])
    ax.set_xlabel("Line number")
    ax.set_ylabel("Time(s)")
    plt.show()


if __name__ == "__main__":
    pts = [(random.random(), random.random()) for k in range(100000)]
    #generate_cutting_size_test(pts, 1000, 5000, 4)
    generate_cutting_time_test(pts, 30000, 40000, 2)
