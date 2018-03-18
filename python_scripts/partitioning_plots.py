import random
import line_testing
import Partitioning
import csv
import matplotlib.pyplot as plt

partitioning_map = {'pmat_poly_pts_color': '#e41a1c',
                    "pmat_poly_lts_color": '#377eb8',
                    "pmat_trap_pts_color": '#f781bf',
                    "pmat_trap_lts_color" : '#4daf4a',
                    "pchan_poly_pts_color": '#984ea3',
                    "pchan_poly_lts_color": '#ff7f00',
                    "pchan_trap_pts_color": '#ffff33',
                    "pchan_trap_lts_color": '#a65628',
                     "sample_color" : '#999999',
                    'pmat_poly_pts_m': 'x',
                    "pmat_poly_lts_m": 'v',
                    "pmat_trap_pts_m": 'o',
                    "pmat_trap_lts_m": '^',
                    "pchan_poly_pts_m": '*',
                    "pchan_poly_lts_m": 'H',
                    "pchan_trap_pts_m": 'D',
                    "pchan_trap_lts_m": '>',
                    "sample_m": '<',
                    'pmat_poly_pts_name': 'Mat_Poly_Points',
                    "pmat_poly_lts_name": 'Mat_Poly_Lines',
                    "pmat_trap_pts_name": 'Mat_Trap_Points',
                    "pmat_trap_lts_name": 'Mat_Trap_Lines',
                    "pchan_poly_pts_name": 'Chan_Poly_Points',
                    "pchan_poly_lts_name": 'Chan_Poly_Lines',
                    "pchan_trap_pts_name": 'Chan_Trap_Points',
                    "pchan_trap_lts_name": 'Chan_Trap_Lines',
                    "sample_name": "Random Sample",
                    "output_size": "Output Size",
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

def plot_partitioning_test_axis(ax, fnames, result_name):

    handles = []
    vparam_name = None
    for fname in fnames:
        with open(fname) as file:
            vparam = []
            results = []


            reader_obj = csv.DictReader(file)
            row = None
            for row in reader_obj:
                if vparam_name != row["vparam"] and vparam_name is not None:
                    raise ValueError("Bad vparam %s in file %s"%(row["vparam"], fname))
                vparam_name = row["vparam"]

                vparam.append(float(row[row["vparam"]]))
                results.append(float(row[result_name]))
            if row is not None:
                if row["part_f"] != "sample":
                    name = row["part_f"] + "_" + row["cutting_f"] + "_" + row["test_set_f"]
                else:
                    name = "sample"
                print(name)
                ax.plot(vparam, results,
                            marker=partitioning_map[name + "_m"],
                            color=partitioning_map[name + "_color"],
                            label=partitioning_map[name + "_name"],
                            linestyle='None')


    #ax.legend(handles=handles)
    ax.legend()
    ax.set_xlabel(partitioning_map[vparam_name])
    ax.set_ylabel(partitioning_map[result_name])


def plot_partitioning(fnames, result_name, output_name):
    ax = plt.subplot()
    plot_partitioning_test_axis(ax, fnames, result_name)
    plt.savefig(output_name + ".png",
                bbox_inches='tight',
                dpi=partitioning_map["dpi"],
                figsize=partitioning_map["shape"]
                )
    plt.show()


fieldnames = ["vparam", "r", "input_size", "output_size", "cell_size",
                  "test_set_f", "cutting_f", "part_f", "time", "error", "k"]

if __name__ == "__main__":
    plot_partitioning(["timing_plot_r2_poly_pts.csv",
                        "timing_plot_r2_poly_lts.csv",
                       "timing_plot_r2_trap_lts.csv",
                        "timing_plot_r2_trap_pts.csv",
                       "sampling.csv"],
                      "time", "time_plot_r2"
                      )

    # plot_partitioning(["timing_plot_r4_poly_pts.csv",
    #                     "timing_plot_r4_poly_pts.csv",
    #                    "timing_plot_r4_trap_lts.csv",
    #                     "timing_plot_r4_trap_pts.csv"],
    #                   "error"
    #                   )



