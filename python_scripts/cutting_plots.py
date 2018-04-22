import random
import numpy as np
import matplotlib.pyplot as plt
import test_set
import poly_tree
import seidel_tree
import poly_tree as pt
import seidel_tree as st
import time
import test_set
import csv

cutting_map = {"poly_d" : "#1b9e77",
               "poly_l" : "#d95f02",
               "poly_p" : "#7570b3",
               "trap_d":"#e7298a",
               "trap_l":"#66a61e",
               "trap_p":"#e6ab02",
               "poly_d_m":"H",
                "poly_l_m":"v",
               "poly_p_m":"o",
               "trap_d_m":">",
               "trap_p_m":"<",
               "trap_l_m": "D"}


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

    ax.plot(line_count, trap_count, marker=marker, c=color)


def plot_cutting_size_b(ax, pts, l_s, h_s, k, test_set_f, b=4, marker='.', color='r'):

    line_count = []
    trap_count = []
    for i in np.linspace(l_s, h_s, k):

        test_set = test_set_f(pts, int(i))
        weight_map = {t: 1 for t in test_set}
        tree = pt.PolyTree2([], test_set)
        tree.cutting_b(b, {l:1 for l in test_set})
        r = len(test_set) / max(cell.total_weight({l:1 for l in test_set}) for cell in tree.get_leaves())
        trap_count.append(b / r**2)
        line_count.append(len(test_set))

    ax.plot(line_count, trap_count, marker=marker, c=color)


def plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f, test_set_f, r=4, marker='.', color='r'):

    line_count = []
    times = []
    for i in np.linspace(l_s, h_s, k):
        print(i)
        t = time.time()
        test_set = test_set_f(pts, int(i))
        weight_map = {t: 1 for t in test_set}

        tree = cutting_f(test_set, weight_map, pts, r)
        e = time.time()
        times.append(e - t)
        line_count.append(len(test_set))

    ax.plot(line_count, times, marker=marker, c=color)


def plot_cutting_size_r(ax, pts, l_s, h_s, k, cutting_f, test_set_f, marker='.', color='r'):

    line_count = []
    trap_count = []
    for i in np.linspace(l_s, h_s, k):
        print(i)
        test_set = test_set_f(pts, int(1/2.0 * len(pts)))
        weight_map = {t: 1 for t in test_set}
        tree = cutting_f(test_set, weight_map, pts, i)
        trap_count.append(len(tree.get_leaves()) / float(i * i))
        line_count.append(i)

    ax.plot(line_count, trap_count, marker=marker, c=color)


def plot_cutting_time_r(ax, pts, l_s, h_s, k, cutting_f, test_set_f, marker='.', color='r'):

    line_count = []
    times = []
    for i in np.linspace(l_s, h_s, k):
        print(i)
        test_set = test_set_f(pts, int(1/2.0 * len(pts)))
        weight_map = {t: 1 for t in test_set}
        t = time.time()
        tree = cutting_f(test_set, weight_map, pts, i)
        e = time.time()
        times.append(e - t)
        line_count.append(i)

    ax.plot(line_count, times, marker=marker, c=color)


def generate_cutting_size_test(pts, l_s, h_s, k, r=4):
    f, ax = plt.subplots()

    print("dual")
    plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_dual, r=r,
                      marker=cutting_map["poly_d_m"],
                      color=cutting_map["poly_d"])
    print("points ")
    plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_points, r=r,
                      marker=cutting_map["poly_p_m"],
                      color=cutting_map["poly_p"])
    print("lines")
    plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_lines, r=r,
                      marker=cutting_map["poly_l_m"],
                      color=cutting_map["poly_l"])

    print("dual")
    plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_dual, r=r,
                      marker=cutting_map["trap_d_m"],
                      color=cutting_map["trap_d"])
    print("points")
    plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_points, r=r,
                      marker=cutting_map["trap_p_m"],
                      color=cutting_map["trap_p"])
    print("lines")
    plot_cutting_size(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_lines, r=r,
                      marker=cutting_map["trap_l_m"],
                      color=cutting_map["trap_l"])
    #ax.legend(["poly dts", "poly psd", "poly lsd", "trap dts", "trap psd", "trap lsd"])
    ax.legend(["PolyTree_8 Dual", "PolyTree_8 Points", "PolyTree_8 Lines",  "Trapezoid Dual", "Trapezoid Points", "Trapezoid Lines"])
    ax.set_xlabel("Number of Lines")
    ax.set_ylabel("Cutting Constant")


    f.savefig("size_test.pdf",
          bbox_inches='tight')


def generate_cutting_size_test_b(pts, l_s, h_s, k, b=4):
    f, ax = plt.subplots()

    print("dual")
    plot_cutting_size_b(ax, pts, l_s, h_s, k, test_set_f=test_set.test_set_dual,  b=b,
                      marker=cutting_map["poly_d_m"],
                      color=cutting_map["poly_d"])
    print("points ")
    plot_cutting_size_b(ax, pts, l_s, h_s, k,  test_set_f=test_set.test_set_points, b=b,
                      marker=cutting_map["poly_p_m"],
                      color=cutting_map["poly_p"])
    print("lines")
    plot_cutting_size_b(ax, pts, l_s, h_s, k, test_set_f=test_set.test_set_lines, b=b,
                      marker=cutting_map["poly_l_m"],
                      color=cutting_map["poly_l"])

    #ax.legend(["poly dts", "poly psd", "poly lsd", "trap dts", "trap psd", "trap lsd"])
    ax.legend(["PolyTree_8 Dual", "PolyTree_8 Points", "PolyTree_8 Lines",  "Trapezoid Dual", "Trapezoid Points", "Trapezoid Lines"])
    ax.set_xlabel("Number of Lines")
    ax.set_ylabel("Cutting Constant")
    plt.show()


def generate_cutting_time_test(pts, l_s, h_s, k, r=4):
    f, ax = plt.subplots()

    print("dual")
    plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_dual, r=r,
                      marker=cutting_map["poly_d_m"],
                      color=cutting_map["poly_d"])
    print("points ")
    plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_points, r=r,
                      marker=cutting_map["poly_p_m"],
                      color=cutting_map["poly_p"])
    print("lines")
    plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_lines, r=r,
                      marker=cutting_map["poly_l_m"],
                      color=cutting_map["poly_l"])

    print("dual")
    plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_dual, r=r,
                      marker=cutting_map["trap_d_m"],
                      color=cutting_map["trap_d"])
    print("points")
    plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_points, r=r,
                      marker=cutting_map["trap_p_m"],
                      color=cutting_map["trap_p"])
    print("lines")
    plot_cutting_time(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_lines, r=r,
                      marker=cutting_map["trap_l_m"],
                      color=cutting_map["trap_l"])
    #ax.legend(["poly dts", "poly psd", "poly lsd", "trap dts", "trap psd", "trap lsd"])
    #ax.legend(["PolyTree_8 Points", "PolyTree_8 Lines",  "Trapezoid Points", "Trapezoid Lines"])
    ax.legend(["PolyTree_8 Dual", "PolyTree_8 Points", "PolyTree_8 Lines",  "Trapezoid Dual", "Trapezoid Points", "Trapezoid Lines"])

    ax.set_xlabel("Number of Lines")
    ax.set_ylabel("Time (sec)")
    f.savefig("time_test.pdf",
                bbox_inches='tight')


def generate_cutting_size_test_r(pts, l_s, h_s, k):
    f, ax = plt.subplots()

    print("dual")
    plot_cutting_size_r(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_dual,
                      marker=cutting_map["poly_d_m"],
                      color=cutting_map["poly_d"])
    print("points ")
    plot_cutting_size_r(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_points,
                      marker=cutting_map["poly_p_m"],
                      color=cutting_map["poly_p"])
    print("lines")
    plot_cutting_size_r(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_lines,
                      marker=cutting_map["poly_l_m"],
                      color=cutting_map["poly_l"])

    print("dual")
    plot_cutting_size_r(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_dual,
                      marker=cutting_map["trap_d_m"],
                      color=cutting_map["trap_d"])
    print("points")
    plot_cutting_size_r(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_points,
                      marker=cutting_map["trap_p_m"],
                      color=cutting_map["trap_p"])
    print("lines")
    plot_cutting_size_r(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_lines,
                      marker=cutting_map["trap_l_m"],
                      color=cutting_map["trap_l"])
    #ax.legend(["poly dts", "poly psd", "poly lsd", "trap dts", "trap psd", "trap lsd"])
    #ax.legend(["PolyTree_8 Points", "PolyTree_8 Lines",  "Trapezoid Points", "Trapezoid Lines"])

    ax.legend(["PolyTree_8 Dual", "PolyTree_8 Points", "PolyTree_8 Lines",  "Trapezoid Dual", "Trapezoid Points", "Trapezoid Lines"])
    ax.set_xlabel("r")
    ax.set_ylabel("Cutting Constant")
    f.savefig("size_test_r.pdf",
                bbox_inches='tight')


def generate_cutting_time_test_r(pts, l_s, h_s, k):
    f, ax = plt.subplots()

    print("dual")
    plot_cutting_time_r(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_dual,
                      marker=cutting_map["poly_d_m"],
                      color=cutting_map["poly_d"])
    print("points ")
    plot_cutting_time_r(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_points,
                      marker=cutting_map["poly_p_m"],
                      color=cutting_map["poly_p"])
    print("lines")
    plot_cutting_time_r(ax, pts, l_s, h_s, k, cutting_f=pt.compute_cutting, test_set_f=test_set.test_set_lines,
                      marker=cutting_map["poly_l_m"],
                      color=cutting_map["poly_l"])

    print("dual")
    plot_cutting_time_r(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_dual,
                      marker=cutting_map["trap_d_m"],
                      color=cutting_map["trap_d"])
    print("points")
    plot_cutting_time_r(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_points,
                      marker=cutting_map["trap_p_m"],
                      color=cutting_map["trap_p"])
    print("lines")
    plot_cutting_time_r(ax, pts, l_s, h_s, k, cutting_f=st.compute_cutting, test_set_f=test_set.test_set_lines,
                      marker=cutting_map["trap_l_m"],
                      color=cutting_map["trap_l"])
    #ax.legend(["poly dts", "poly psd", "poly lsd", "trap dts", "trap psd", "trap lsd"])
    #ax.legend(["PolyTree_8 Points", "PolyTree_8 Lines",  "Trapezoid Points", "Trapezoid Lines"])
    ax.legend(["PolyTree_8 Dual", "PolyTree_8 Points", "PolyTree_8 Lines",  "Trapezoid Dual", "Trapezoid Points", "Trapezoid Lines"], loc="upper left")

    ax.set_xlabel("r")
    ax.set_ylabel("Time (sec)")
    f.savefig("time_test_r.pdf",
                bbox_inches='tight')


if __name__ == "__main__":
    pts = [(random.random(), random.random()) for k in range(10000)]
    # generate_cutting_size_test_r(pts, 2, 32, 20)
    #generate_cutting_time_test_r(pts, 2, 32, 20)


    #generate_cutting_time_test_b(pts, 2, 30, 10)
    generate_cutting_size_test(pts, 200, 2000, 20)
    generate_cutting_time_test(pts, 200, 2000, 20)
