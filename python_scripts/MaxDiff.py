import pyscan
import argparse
import random
import csv
import bisect
import time
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import Counter

def sumMatchingR(rows):
    r_vals = set((row['r'] for row in rows))
    for r in r_vals:
        alg_time_count = {'ff': {}, 'fs': {}, 'ss':{}, 'ns':{}}
        r_counts = {'ff': {}, 'fs': {}, 'ss':{}, 'ns':{}}
        for row in rows:
            if row['r'] in alg_time_count[row['alg']]:
                alg_time_count[row['alg']][row['r']][0] += float(row['maxDiff'])
                alg_time_count[row['alg']][row['r']][1] += float(row['time'])
                r_counts[row['alg']][row['r']] += 1
            else:
                alg_time_count[row['alg']][row['r']] = [float(row['maxDiff']), float(row['time'])]
                r_counts[row['alg']][row['r']] = 1

        for alg in r_counts:
            for r in r_counts[alg]:
                alg_time_count[alg][r][0] = alg_time_count[alg][r][0] / r_counts[alg][r]
                alg_time_count[alg][r][1] = float(r) #alg_time_count[alg][r][1] / r_counts[alg][r]
        return alg_time_count

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--prefix',
                        required=True,
                        help="prefix for the timing and error outputs")
    args = parser.parse_args()
    with open(args.prefix + "_etime.csv", 'r') as csvFile:
        reader = csv.DictReader(csvFile)
        all_rows = [row for row in reader]
        row_sums = sumMatchingR(all_rows)
        #print row_sums

        time_ff = [float(row['time']) for row in all_rows if row['alg'] == 'ff']
        time_fs = [float(row['time']) for row in all_rows if row['alg'] == 'fs']
        time_ss = [float(row['time']) for row in all_rows if row['alg'] == 'ss']
        time_ns = [float(row['time']) for row in all_rows if row['alg'] == 'ns']

        err_ff = [float(row['maxDiff']) for row in all_rows if row['alg'] == 'ff']
        err_fs = [float(row['maxDiff']) for row in all_rows if row['alg'] == 'fs']
        err_ss = [float(row['maxDiff']) for row in all_rows if row['alg'] == 'ss']
        err_ns = [float(row['maxDiff']) for row in all_rows if row['alg'] == 'ns']

        """
        ff = [row_sums['ff'][el] for el in row_sums['ff']]
        fs = [row_sums['fs'][el] for el in row_sums['fs']]
        ss = [row_sums['ss'][el] for el in row_sums['ss']]
        ns = [row_sums['ns'][el] for el in row_sums['ns']]

        ff = sorted(ff, key=lambda x: x[1])
        fs = sorted(fs, key=lambda x: x[1])
        ss = sorted(ss, key=lambda x: x[1])
        ns = sorted(ns, key=lambda x: x[1])
        err_ff, time_ff = zip(*ff)
        #err_fs, time_fs = zip(*fs)
        err_ss, time_ss = zip(*ss)
        err_ns, time_ns = zip(*ns)
        """
        #print row_sums['ff']
        print err_ff
        print time_ff
        plt.hold(True)
        ax = plt.subplot(1,1,1)
        p1, = ax.plot(time_ff, err_ff, 'o', color='g', label='SlabScan')
        p2, = ax.plot(time_fs, err_fs, 'o', color='b', label='SubSumScan')

        if err_ss != []:
            p3, = ax.plot(time_ss, err_ss, '.', color='r', label='GridScan')
            p4, = ax.plot(time_ns, err_ns, 'x', color='k', label='NetScan')


        ax.legend(loc='upper right')
        ax.set_xlabel('Time (sec)')
        ax.set_ylabel('Diff from Planted Max')
        plt.show()