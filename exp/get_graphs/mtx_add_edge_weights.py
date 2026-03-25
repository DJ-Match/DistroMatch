import argparse
import random
import numpy as np
import tempfile
import shutil
import os
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter

MAX_WEIGHT = 500_000
MIN_WEIGHT = 1
NUM_THREADS = 8  # You can tweak based on your CPU

def generate_uniform_weight():
    return random.randint(MIN_WEIGHT, MAX_WEIGHT)

def generate_exponential_weight(scale=100000.0):
    w = np.random.exponential(scale)
    return int(min(MAX_WEIGHT, max(MIN_WEIGHT, w)))

def split_header_and_data(lines):
    header_lines = []
    size_line = None
    data_lines = []

    for line in lines:
        if line.startswith('%'):
            # Replace pattern line with integer
            if line.strip().lower() == '%%matrixmarket matrix coordinate pattern symmetric':
                line = '%%MatrixMarket matrix coordinate integer symmetric\n'
            header_lines.append(line)
        elif size_line is None and line.strip():
            size_line = line
        else:
            data_lines.append(line)

    return header_lines, size_line, data_lines

def process_file_streaming(input_file, output_file, weight_fn, plot=False, scale=100000.0):
    weights_counter = Counter() if plot else None

    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        header_lines = []
        size_line = None

        # First pass: handle headers
        for line in fin:
            if line.startswith('%'):
                if line.strip().lower() == '%%matrixmarket matrix coordinate pattern symmetric':
                    line = '%%MatrixMarket matrix coordinate integer symmetric\n'
                header_lines.append(line)
            else:
                size_line = line
                break

        # Write header
        fout.writelines(header_lines)
        if size_line:
            fout.write(size_line)

        # Stream rest of the file
        for line in fin:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            u, v = parts[:2]
            w = weight_fn() if weight_fn != generate_exponential_weight else generate_exponential_weight(scale)
            fout.write(f"{u} {v} {w}\n")
            if plot:
                weights_counter[w] += 1

    return weights_counter



def plot_weight_frequency(counter, title, filename):
    x = list(counter.keys())
    y = list(counter.values())

    plt.figure(figsize=(12, 6))
    plt.scatter(x, y, s=12, alpha=0.7, color='mediumseagreen')
    plt.title(title)
    plt.xlabel("Edge Weight")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Add random integer edge weights to MTX graph.")
    parser.add_argument("input", help="Input MTX file")
    parser.add_argument("--exp_scale", type=float, default=30000.0, help="Scale for exponential distribution")
    parser.add_argument("--plot", action="store_true", help="Plot histogram of edge weights")
    args = parser.parse_args()

    input_path = args.input
    base, ext = os.path.splitext(input_path)
    uniform_out = f"{base}_uni{ext}"
    exp_out = f"{base}_exp{ext}"

    print(f"Processing uniform weights → {uniform_out}")
    uniform_counter = process_file_streaming(input_path, uniform_out, generate_uniform_weight, plot=args.plot)

    print(f"Processing exponential weights in-place")
    fd, tmp_path = tempfile.mkstemp()
    os.close(fd)  # Close immediately, will reopen
    exp_counter = process_file_streaming(input_path, tmp_path, generate_exponential_weight, plot=args.plot, scale=args.exp_scale)


    print(f"Renaming to → {exp_out}")
    shutil.move(tmp_path, exp_out)

    if args.plot:
        print("Generating plots...")
        plot_weight_frequency(uniform_counter, "Uniform Edge Weight Frequency", f"{base}_uni_scatter.png")
        plot_weight_frequency(exp_counter, "Exponential Edge Weight Frequency", f"{base}_exp_scatter.png")



    print("Done.")

if __name__ == "__main__":
    main()
