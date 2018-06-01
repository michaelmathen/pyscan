import argparse
import csv
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--prefix',
                        required=True,
                        help="prefix for the timing and error outputs")
    args = parser.parse_args()
    with open(args.prefix + "_error.csv", 'r') as csvFile:
        reader = csv.DictReader(csvFile)
        all_rows = [row for row in reader]
        error = [row['error'] for row in all_rows]

        r_prime = [row['r_prime'] for row in all_rows]
        plt.hold(True)
        ax = plt.subplot(1,1,1)
        p, = ax.plot(r_prime, error, label='Function Error')
        ax.legend(loc='upper right')
        ax.set_xlabel('r\'')
        ax.set_ylabel('Error')
        plt.show()