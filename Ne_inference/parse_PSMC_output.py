#!/usr/bin/env python3

import os
import re
import numpy as np
import matplotlib.pyplot as plt

def weighted_mean(data, weights):
    return np.sum(np.array(data) * np.array(weights)) / np.sum(weights)

def weighted_geometric_mean(data, weights):
    log_data = np.log(data)
    return np.exp(np.sum(log_data * weights) / np.sum(weights))

def weighted_harmonic_mean(data, weights):
    return np.sum(weights) / np.sum(weights / np.array(data))

class PSMCPlotter:
    def __init__(self, out_prefix, psmc_file, mutation_rate=2.5e-8, skip=100, max_gen=0, min_gen=10000,
                 max_popsize=0, min_iterations=5, nth_iteration=20, multiline_titles='', false_neg_rate=0.0,
                 years_per_gen=25, font='Helvetica,16', line_width=4, key_position='right top', title='',
                 no_scaling=False, show_last_bin=False, convert_to_pdf=False, keep_files=False, plot_grid=False):
        self.out_prefix = out_prefix
        self.psmc_file = psmc_file
        self.mutation_rate = mutation_rate
        self.skip = skip
        self.max_gen = max_gen
        self.min_gen = min_gen
        self.max_popsize = max_popsize
        self.min_iterations = min_iterations
        self.nth_iteration = nth_iteration
        self.multiline_titles = multiline_titles
        self.false_neg_rate = false_neg_rate
        self.years_per_gen = years_per_gen
        self.font = font
        self.line_width = line_width
        self.key_position = key_position
        self.title = title
        self.no_scaling = no_scaling
        self.show_last_bin = show_last_bin
        self.convert_to_pdf = convert_to_pdf
        self.keep_files = keep_files
        self.plot_grid = plot_grid
        self.data = []
        self.modifiers = {}

    def parse_pattern(self, pattern):
        terms = pattern.split('+')
        n_lambda = 0
        stack = []
        for term in terms:
            parts = term.split('*')
            count, value = (int(parts[0]), int(parts[1])) if len(parts) == 2 else (1, int(parts[0]))
            stack.extend([value] * count)
            n_lambda += count
        ret = [0] * sum(stack)
        pos = 0
        for num in stack:
            for _ in range(num):
                ret[pos] = num
                pos += 1
        return n_lambda, ret

    def load_data(self):
        max_intv = 10000
        id = 0
        do_store = False
        gof = 'RI'
        min_ri = float('inf')
        N_scale = 1.0
        Mseg = Msize = dt = 0
        data = []
        modifiers = {'FN': [], 'nscale': [], 'tscale': [], 'alpha': [], 'xshift': []}
        with open(self.psmc_file, 'r') as file:
            for line in file:
                if match := re.search(r'^MM.*skip:(\d+)', line):
                    self.skip = int(match.group(1)) * self.skip
                elif match := re.search(r'^MM.*pattern:(\S+),', line):
                    n_lambda, pat = self.parse_pattern(match.group(1))
                    for i in reversed(range(len(pat))):
                        if pat[i] != pat[-1]: break
                    max_intv = i
                elif re.search(r'^MM.*is_decoding:', line):
                    d = {'T': 0, 'R': 0, 'N0': 0, 'RI': 0, 'D': []}
                    data.append(d)
                    id += 1
                    min_ri = float('inf')
                elif match := re.search(r'^RD\s(\S+)', line):
                    round = int(match.group(1))
                elif (match := re.search(r'^(RI|GF)\s(\S+)', line)) and match.group(1) == gof and 'nan' not in match.group(2) and 'inf' not in match.group(2):
                    value = float(match.group(2))
                    do_store = round >= self.min_iterations and value < min_ri
                    if self.nth_iteration > 0:
                        do_store = (round == self.nth_iteration)
                elif do_store and (match := re.search(r'^TR\s(\S+)\s(\S+)', line)):
                    d['T'], d['R'] = float(match.group(1)) / self.skip, float(match.group(2)) / self.skip
                    d['N0'] = d['T'] / (4 * self.mutation_rate) * N_scale
                    d['RI'] = min_ri
                elif do_store and (match := re.search(r'^DT\s(\S+)', line)):
                    dt = float(match.group(1))
                elif do_store and (match := re.search(r'^RS\s(\d+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)', line)):
                    if self.show_last_bin or int(match.group(1)) <= max_intv:
                        seg, scaleN0, gen4, gen5, gen6 = map(float, match.groups()[1:])
                        d['D'].append((2 * d['N0'] * (seg + dt) * self.years_per_gen, scaleN0 * d['N0'] / 10000, gen4, gen5, gen6))
                        Mseg = max(Mseg, gen4)
                        Msize = max(Msize, 2 * d['N0'] * seg)
                elif do_store and (match := re.search(r'^PA\s(.*)', line)):
                    d['PAR'] = match.group(1)
                elif do_store and re.search(r'^//', line):
                    d['Mseg'] = Mseg
                    d['Msize'] = Msize
        self.data = data

    def output_stats(self, label, statsfile):
        self.load_data()
        t_values = []
        n_values = []
        for d in self.data:
            t_values = [entry[0] for entry in d['D']]
            n_values = [entry[1] for entry in d['D']]
        if not self.no_scaling:
            k = next((i for i, x in enumerate(t_values) if x > self.min_gen), 0)
            t_values = t_values[k:]
            n_values = n_values[k:]
        mean = sum(n_values) / len(n_values)
        geom_mean = np.exp(sum(np.log(n_values)) / len(n_values))
        harm_mean = len(n_values) / sum(1/y for y in n_values)
        weights = np.diff(t_values, prepend=0)
        w_mean = weighted_mean(n_values, weights)
        w_geom_mean = weighted_geometric_mean(n_values, weights)
        w_harm_mean = weighted_harmonic_mean(n_values, weights)
        with open(statsfile, 'a') as file:
            file.write(f"{label}\t{mean}\t{geom_mean}\t{harm_mean}\t{w_mean}\t{w_geom_mean}\t{w_harm_mean}\n")
        return t_values, n_values

if __name__ == "__main__":
    os.makedirs("output", exist_ok=True)
    with open("output/no_scale_stats.tsv", "w") as file:
        file.write("label\tmean\tgeom_mean\tharm_mean\tw_mean\tw_geom_mean\tw_harm_mean\n")
    with open("output/scale_stats.tsv", "w") as file:
        file.write("label\tmean\tgeom_mean\tharm_mean\tw_mean\tw_geom_mean\tw_harm_mean\n")

    for file_name in os.listdir("illumina_sim/psmc_out"):
        if file_name.startswith("."): continue
        file_path = os.path.join("illumina_sim/psmc_out", file_name)
        label = file_name.rsplit(".", 1)[0]

        plotter1 = PSMCPlotter("example_output", file_path, no_scaling=True)
        plotter2 = PSMCPlotter("example_output", file_path, no_scaling=False)

        plotter1.output_stats(label, "output/no_scale_stats.tsv")
        plotter2.output_stats(label, "output/scale_stats.tsv")
