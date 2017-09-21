import pyscan
import argparse
import random
import csv
import bisect
import time
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import Counter

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--prefix',
                        required=True,
                        help="prefix for the timing and error outputs")
    args = parser.parse_args()
    with open(args.prefix + "_time.csv", 'r') as csvFile:
        reader = csv.DictReader(csvFile)
        all_rows = [row for row in reader]
        time_ff = [row['time'] for row in all_rows if row['alg'] == 'ff']
        time_fs = [row['time'] for row in all_rows if row['alg'] == 'fs']
        time_ss = [row['time'] for row in all_rows if row['alg'] == 'ss']

        r_ff = [row['r'] for row in all_rows if row['alg'] == 'ff']
        r_fs = [row['r'] for row in all_rows if row['alg'] == 'fs']
        r_ss = [row['r'] for row in all_rows if row['alg'] == 'ss']
        plt.hold(True)
        ax = plt.subplot(1,1,1)
        p1, = ax.plot(r_ff, time_ff, color='g', label='SlabScan')
        #p2, = ax.plot(r_fs, time_fs, color='b', label='SubSumScan')
        if r_ss != []:
            p3, = ax.plot(r_ss, time_ss, color='r', label='GridScan/NetScan')
        ax.legend(loc='upper left')
        ax.set_xlabel('Grid Size')
        ax.set_ylabel('Time (sec)')
        plt.show()
