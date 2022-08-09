import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import fnmatch
from pathlib import Path
from Bio import SeqIO
import heapq
import time

from parameters import info_file, fig_size_inches, genome_bin_size, low_reads_cap_percent, plots_filetype, plots_dpi
from .utils import output_path, get_log_entry

####
# Makes 3 graphs:
# 1. A binned histogram of reads
# 2. A modified version of the first binned histogram with y-axis capped at a custom percentage, to show missed reads
# 3. A zoomed-in histogram around the target region, showing nearing misses
###

# font control
plt.rcParams['svg.fonttype'] = 'none'  # important so that text stays as intact characters in the output
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

plt.ioff()

def get_bins(filename, genome_length, bin_size):
    csv = np.genfromtxt(filename, delimiter=",", skip_header=1)
    if len(csv.shape) == 1:
        csv = np.array([csv.tolist()])
    number_of_bins = (genome_length // bin_size) + 1
    bin_numbers = np.array(range(1, number_of_bins + 1))
    bin_values = []
    k = 0
    for i in bin_numbers:
        min_val = (i - 1) * bin_size
        max_val = i * bin_size
        bin_counts = 0
        while (k < len(csv) and csv[k][0] <= max_val):
            bin_counts += csv[k][1]
            k += 1
            if k == len(csv):
                break
        #counts_in_bin = [entry[1] for entry in csv if min_val <= entry[0] < max_val]
        #total_for_bin = sum(counts_in_bin)
        bin_values.append(bin_counts)
    np_bin_values = np.asarray(bin_values)
    return bin_numbers, np_bin_values

def setup_axes(axs, max_x, max_y, genome_length, bin_size):
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.spines['bottom'].set_position('zero')
    # configure spines of the graph; 'remove' top and right spines
    
    axs.spines['left'].set_bounds(0, max_y)
    axs.spines['bottom'].set_bounds(0, max_x)  # bottom line starts and end with length of genome
    axs.spines['left'].set_position('zero')

    # axis labeling
    #plt.xlabel("E. coli genomic coordinate (Mbp)")
    #plt.ylabel("Read count")
    axs.yaxis.set_label_coords(-0.02, 0.5)

    # set up xticks
    increments = [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 100000000, 1000000000, 10000000000]
    i = 0
    increment = increments[0]
    while genome_length // increment > 10:
        i += 1
        increment = increments[i]
    
    x_axis_size = int(genome_length / bin_size)

    axs.set_xticks(range(0, x_axis_size, int(increment / bin_size)))
    axs.set_xticklabels(np.arange(0, genome_length / 1000000, increment / 1000000))

    # set limits of graph area - different from spine limits
    # leave extra space on top for info text, and space on bottom for the spacer marker
    plt.gca().set_xlim(left=-0.15*max_x, right=max_x)
    plt.gca().set_ylim(bottom=-0.12*max_y, top=1.20*max_y)

    return axs

def plot_binned(filepath, run_information, yAxis_type, isPlasmid=False):  
    # Main graphing function for histogram binning
    # The yAxis_type can one of either:
    # 1. "raw" - the raw counts of reads in each bin
    # 2. "normalized" - by percentage of total reads
    # 3. "zoomed" - a normalized graph that is zoomed in at the y-axis to show low-frequence bins
    code = run_information.get('Sample', None)
    description = run_information.get('Description', None)
    exp_date = run_information.get('Experiment date', None)
    spacer_locations = [int(loc) for loc in run_information['End of protospacer'].split()]
    genome_path = run_information['Target fasta file']
    if isPlasmid:
        genome_path = run_information['Second target fasta file']
    genome_length = len(SeqIO.read(Path(genome_path), 'fasta'))
    bin_scale = int(100000/genome_bin_size)
    bin_size = genome_bin_size
    if isPlasmid or genome_length < 100000:
        bin_scale = int(genome_length // 200)
        bin_size = 100
    
    # determine which bin the spacer lies in
    spacer_bins = []
    if not isPlasmid:
        spacer_bins = [int(spacer_location / bin_size) + 0.5 for spacer_location in spacer_locations]
    # get bins and counts for the histogram
    b, a2 = get_bins(filepath, genome_length, bin_size)
    total_reads = int(sum(a2))

    max_x = len(a2)

    # For normalized and zoomed graphs, normalized based on total reads as a percentage
    max_y = max_y_label = None
    if yAxis_type == 'zoomed':
        max_y = total_reads * low_reads_cap_percent * 0.01
        max_y_label = low_reads_cap_percent
        a2 = np.fmin(a2, max_y)
    elif yAxis_type == 'normalized':
        max_y = 100
        max_y_label = '100%'
        a2 = (100*a2)/total_reads
    elif yAxis_type == 'raw':
        max_y = int(max(a2))
        max_y_label = str(max_y)
    else:
        print('Must give a yAxis type or raw, normalized, or zoomed')

    # set up figure
    fig, axs = plt.subplots(1, 1, tight_layout=True)
    setup_axes(axs, max_x, max_y, genome_length, bin_size)
    # 2 ticks on the y axis, one at 0 and on at max value of y (max reads)
    axs.set_yticks([0, max_y])
    axs.set_yticklabels([0, max_y_label])
    
    # 1. scatter plot points for protospacer markers
    for spacer_bin in spacer_bins:
        axs.scatter(spacer_bin, -0.03 * max_y, marker="^", c="#A32E79", s=33)
    # 2. main bar graph, width slightly >1 to avoid gaps between bars
    axs.bar(b, a2, color='#0D62AC', width=1.05)  

    # size of output figures
    fig.set_size_inches(fig_size_inches[0], fig_size_inches[1])

    # The text on the graph
    text_x = 5 * bin_scale
    text_y = 1.17 * max_y
    if yAxis_type == 'zoomed':
        text_x = 33 * bin_scale
        text_y = 0.8 * max_y
    plt.text(text_x, text_y, "{}-{}\n{}\nTotal Reads = {}"
        .format(code, exp_date, description, total_reads))

    # save figure
    sample = run_information['Sample']
    graph_name = "plasmid" if isPlasmid else "genome"
    plt.savefig(output_path(os.path.join('plots', '{}_{}_hist_{}.{}'.format(sample, graph_name, yAxis_type, plots_filetype))), dpi=plots_dpi)
    plt.close()  # closes the matplotlib preview popup window

def make_genome_plots(csvFile, meta_info, isPlasmid=False):
    start = time.perf_counter()
    graph_name = "plasmid" if isPlasmid else "genome"
    print("Making {} mapping histograms...".format(graph_name))

    plot_binned(csvFile, meta_info, 'raw', isPlasmid)
    plot_binned(csvFile, meta_info, 'normalized', isPlasmid)
    plot_binned(csvFile, meta_info, 'zoomed', isPlasmid)
    elapsed_time = round(time.perf_counter() - start, 2)
    print("Finished {} mapping plotting in {} seconds".format(graph_name, elapsed_time))
    return meta_info
