import csv
import matplotlib.pyplot as plt
import math
import plotting_tools
import numpy as np
import matplotlib.ticker as tick

partitioning_map = {'pmat_poly_pts_color': '#e41a1c',
                    "pmat_poly_lts_color": '#377eb8',
                    "pmat_poly_dts_color": '#f781bf',
                    "pmat_trap_lts_color" : '#4daf4a',
                    "l2_color": '#984ea3',
                    #"sample_color": '#ff7f00',
                    "pchan_color": '#ff7f00',
                    "pchans_color": '#ffff33',
                    "pchan_trap_lts_color": '#a65628',
                    "box_color": "#6a3d9a",
                     "sample_color" : '#999999',
                    'pmat_poly_pts_m': 'x',
                    "pmat_poly_lts_m": 'v',
                    "pmat_poly_dts_m": 'o',
                    "pchan_m": '.',
                    "pchans_m": 'p',
                    "box_m": ".",
                    "sample_m": '<',
                    "l2_m": ">",
                    'pmat_poly_pts_name': 'Mat Poly Points',
                    "pmat_poly_lts_name": 'Mat Poly Lines',
                    "pmat_poly_dts_name": 'Mat Poly Dual',
                    "l2_name": 'Biased L2',
                    "pmat_trap_lts_name": 'Mat Trap Lines',
                    "pchans_name": 'Chan Simple',
                    "pchan_name": 'Chan',
                    "pchan_poly_lts_name": 'Chan Poly Lines',
                    "pchan_trap_pts_name": 'Chan Trap Points',
                    "pchan_trap_lts_name": 'Chan Trap Lines',
                    "box_name" : "KD Tree Sample",
                    "sample_name": "Random Sample",
                    "output_size": "Output Size",
                    "input_size": "Input Size",
                    "m_disc" : "Discrepancy",
                    "b": "b",
                    "time": "Time (sec)",
                    "error": "Error",
                    "dpi":600,
                    "shape": (8,6)
                    }


def plot_order(names):
    return names.sorted()

# Plots showing sample size vs. error
# Plots showing sample size vs. time
# Plots showing test set size vs. error
# Plots showing r vs. error
# Plots showing test set size vs. time
# Plots showing r vs. time

def plot_partitioning_test_axis(ax, fnames, result_name, vparam_name):
    # if result_name == "error":
    #     ax.set_yscale("log")
        # if vparam_name is "output_size":
        #     ax.set_xscale("log")
    min_result = math.inf
    max_result = -math.inf
    min_vparam = math.inf
    max_vparam = -math.inf

    handles = []
    for fname in fnames:
        with open(fname) as file:
            vparam = []
            results = []


            reader_obj = csv.DictReader(file)
            row = None
            for row in reader_obj:
                # if vparam_name != row["vparam"] and vparam_name is not None:
                #     raise ValueError("Bad vparam %s in file %s"%(row["vparam"], fname))
                # vparam_name = row["vparam"]

                vparam.append(float(row[vparam_name]))
                results.append(float(row[result_name]))

            indices = sorted(range(len(vparam)), key = lambda x: vparam[x])
            vparam = [vparam[i] for i in indices]
            results = [results[i] for i in indices]
            min_vparam = min(min(vparam), min_vparam)
            max_vparam = max(max(vparam), max_vparam)

            min_result = min(min(results), min_result)
            max_result = max(max(results), max_result)
            if row is not None:
                if row["part_f"] == "pmat":
                    name = row["part_f"] + "_" + row["cutting_f"] + "_" + row["test_set_f"]
                else:
                    name = row["part_f"]
                if result_name == "error":
                    #plt.yscale("log")
                    # plt.xscale("log")
                    # ax.yaxis.set_major_locator(tick.LogLocator(base=10.0,numticks=4, subs=(0.001, )))
                    # ax.yaxis.set_major_formatter(tick.LogFormatter(labelOnlyBase=False))
                    ax.plot(vparam, results,
                                marker=partitioning_map[name + "_m"],
                                color=partitioning_map[name + "_color"],
                                label=partitioning_map[name + "_name"])

                else:
                    ax.plot(vparam, results,
                                marker=partitioning_map[name + "_m"],
                                color=partitioning_map[name + "_color"],
                                label=partitioning_map[name + "_name"]
                                )
    ax.set_xlim([min_vparam, max_vparam])
    ax.set_ylim([0, max_result])

    #ax.legend(handles=handles)
    ax.legend(loc="best")
    ax.set_xlabel(partitioning_map[vparam_name])
    ax.set_ylabel(partitioning_map[result_name])


def plot_partitioning(fnames, result_name, output_name, vparam_name, x_min, x_max):
    f, ax = plt.subplots()
    plot_partitioning_test_axis(ax, fnames, result_name, vparam_name)
    # ax.set_xlim([x_min, x_max])
    f.savefig(output_name,
                bbox_inches='tight',
                dpi=partitioning_map["dpi"],
                figsize=partitioning_map["shape"]
                )
    plt.show()



def plot_discrepancy(fnames, names, result_name, vparam_name, output_name, log_x=False, log_y=False):
    f, ax = plt.subplots()

    # if log_y:
    #     ax.set_yscale("log")
    # if log_x:
    #ax.set_xscale("log")

    min_result = math.inf
    max_result = -math.inf
    min_vparam = math.inf
    max_vparam = -math.inf


    for fname, name in zip(fnames, names):
        vparam, results = plotting_tools.read_csv(fname, vparam_name, result_name)
        indices = sorted(range(len(vparam)), key = lambda x: vparam[x])

        vparam = [vparam[i] for i in indices]
        results = [.0134 - results[i] for i in indices]

        plotting_tools.kernel_error_plot(ax, np.log10(vparam), results, vparam_name, result_name,
                                         marker=partitioning_map[name + "_m"],
                                         color=partitioning_map[name + "_color"],
                                         label=partitioning_map[name + "_name"], sig1=.15, alpha=.2, error_bars=False)

        # ax.plot(np.log10(vparam), results,
        #             marker=partitioning_map[name + "_m"],
        #             color=partitioning_map[name + "_color"],
        #             label=partitioning_map[name + "_name"]
        #             )
    #ax.plot([-2, 3], [.0134, .0134], '-', label="Max Discrepancy", linestyle="dashed")
    ax.get_xaxis().set_major_formatter(tick.FormatStrFormatter("$10^{%.1f}$"))
    ax.set_xlabel(partitioning_map[vparam_name])
    ax.set_ylabel("Discrepancy Error")
    ax.legend(loc="best")
    ax.set_xlim([-.2, 1])
    ax.set_ylim([0, .012])
    f.savefig(output_name,
                bbox_inches='tight',
                dpi=partitioning_map["dpi"],
                figsize=partitioning_map["shape"]
                )
    plt.show()

fieldnames = ["vparam", "r", "input_size", "output_size", "cell_size",
                  "test_set_f", "cutting_f", "part_f", "time", "error", "k"]

#
# plot_partitioning(["test.csv",
#                    "output_size_pchan_poly_lts__2_100000_200_1",
#                    "output_sampling_chan.csv"],
#                   "error",
#                   "test.png"
#                   )

if __name__ == "__main__":


    input_directory = "small_comparisons5/"
    # plot_partitioning([input_directory + "input_size_pmat_poly_dts.csv",
    #                    input_directory + "input_size_pmat_poly_lts.csv",
    #                    input_directory + "input_size_pmat_poly_pts.csv",
    #                    input_directory + "input_size_box.csv",
    #                    input_directory + "input_size_pchan_poly_dts.csv"],
    #                   "error",
    #                   input_directory + "input_size_error_plots.pdf"
    #                   )

    # for i in ["input_size"]:
    #     #i = "output_size"
    #     for j in ["error", "time"]:
    #         plot_partitioning([input_directory + i +"_pmat_poly_dts.csv",
    #                            input_directory + i + "_pmat_poly_lts.csv",
    #                            input_directory + i + "_pmat_poly_pts.csv",
    #                            input_directory + i + "_box.csv",
    #                            input_directory + i + "_pchan_poly_dts.csv",
    #                             input_directory + i + "_pchans_poly_dts.csv",
    #                            input_directory + i + "_sample.csv"
    #                             ],
    #                           j,
    #                           input_directory + i + "_" + j + "_plots.pdf",
    #                           i
    #                           ,16,
    #                           128)


    input_directory = "Discrepancy5/"
    plot_discrepancy([
              input_directory + "chan_discrepancy.csv",
               input_directory + "chans_discrepancy.csv",
               input_directory + "mat_discrepancy.csv",
               input_directory + "naive_discrepancy.csv",
               input_directory + "quad_discrepancy.csv"
               ],
              ["pchan",
               "pchans",
               "pmat_poly_pts",
               "sample",
               "box"],
                "m_disc",
              "time",
                input_directory + "s" + "_plots.pdf",
              log_y=False,
              log_x=True)


    # plot_partitioning(["input_size_pchan_poly_lts__2_100000_200_1",
    #                    "input_size_pchan_poly_pts__2_100000_200_1",
    #                     "input_size_pchan_trap_lts__2_100000_200_1",
    #                     "input_size_pchan_trap_pts__2_100000_200_1",
    #                    "input_size_box_100000_200_1",
    #                    "input_sampling_chan.csv"],
    #                   "time",
    #                   "input_size_time_pchan_plot_2_100000_200_1"
    #                   )
    #
    # plot_partitioning(["input_size_pchan_poly_lts__2_100000_200_1",
    #                    "input_size_pchan_poly_pts__2_100000_200_1",
    #                     "input_size_pchan_trap_lts__2_100000_200_1",
    #                     "input_size_pchan_trap_pts__2_100000_200_1",
    #                    "input_size_box_100000_200_1",
    #                    "input_sampling_chan.csv"],
    #                   "error",
    #                   "input_size_error_pchan_plot_2_100000_200_1"
    #                   )

    # plot_partitioning(["r_pchan_poly_lts__100000_200_1",
    #                    "r_pchan_poly_pts__100000_200_1",
    #                    "r_pchan_trap_lts__100000_200_1",
    #                    "r_pchan_trap_pts__100000_200_1"],
    #                   "time",
    #                   "r_time_pchan_plot_100000_200_1"
    #                   )
    #
    # plot_partitioning(["r_pchan_poly_lts__100000_200_1",
    #                    "r_pchan_poly_pts__100000_200_1",
    #                    "r_pchan_trap_lts__100000_200_1",
    #                    "r_pchan_trap_pts__100000_200_1"],
    #                   "error",
    #                   "r_error_pchan_plot_2_100000_200_1"
    #                   )
    #
    #
    #
    # plot_partitioning(["output_size_pmat_poly_lts__4_10000_200_1",
    #                    "output_size_pmat_poly_pts__4_10000_200_1",
    #                    "output_size_pmat_trap_lts__4_10000_200_1",
    #                    "output_size_pmat_trap_pts__4_10000_200_1",
    #                    "output_sampling_mat.csv"],
    #                   "time",
    #                   "output_size_time_pmat_plot_4_10000_200_1"
    #                   )
    #
    # plot_partitioning(["output_size_pmat_poly_lts__4_10000_200_1",
    #                    "output_size_pmat_poly_pts__4_10000_200_1",
    #                    "output_size_pmat_trap_lts__4_10000_200_1",
    #                    "output_size_pmat_trap_pts__4_10000_200_1",
    #                    "output_sampling_mat.csv"],
    #                   "error",
    #                   "output_size_error_pmat_plot_4_10000_200_1"
    #                   )
    #
    # plot_partitioning(["input_size_pmat_poly_lts__4_10000_200_1",
    #                    "input_size_pmat_poly_pts__4_10000_200_1",
    #                    "input_size_pmat_trap_lts__4_10000_200_1",
    #                    "input_size_pmat_trap_pts__4_10000_200_1",
    #                    "input_sampling_mat.csv"],
    #                   "time",
    #                   "input_size_time_pmat_plot_4_10000_200_1"
    #                   )
    #
    # plot_partitioning(["input_size_pmat_poly_lts__4_10000_200_1",
    #                    "input_size_pmat_poly_pts__4_10000_200_1",
    #                    "input_size_pmat_trap_lts__4_10000_200_1",
    #                    "input_size_pmat_trap_pts__4_10000_200_1",
    #                    "input_sampling_mat.csv"],
    #                   "error",
    #                   "input_size_error_pmat_plot_2_10000_200_1"
    #                   )
    #
    # plot_partitioning(["r_pmat_poly_lts__10000_200_1",
    #                    "r_pmat_poly_pts__10000_200_1",
    #                    "r_pmat_trap_lts__10000_200_1",
    #                    "r_pmat_trap_pts__10000_200_1"],
    #                   "time",
    #                   "r_time_pmat_plot_10000_200_1"
    #                   )
    #
    # plot_partitioning(["r_pmat_poly_lts__10000_200_1",
    #                    "r_pmat_poly_pts__10000_200_1",
    #                    "r_pmat_trap_lts__10000_200_1",
    #                    "r_pmat_trap_pts__10000_200_1"],
    #                   "error",
    #                   "r_error_pmat_plot_2_10000_200_1"
    #                   )

#     # plot_partitioning(["timing_plot_r4_poly_pts.csv",
#     #                     "timing_plot_r4_poly_pts.csv",
#     #                    "timing_plot_r4_trap_lts.csv",
#     #                     "timing_plot_r4_trap_pts.csv"],
#     #                   "error"
#     #                   )



