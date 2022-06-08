from __future__ import absolute_import
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
from pipeline.utils import output_path, get_log_entry

plt.rcParams['svg.fonttype'] = 'none'  # important so that text stays as intact characters in the output
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

plt.ioff()

def setup_axes(axs, min_x, max_x, max_y, genome_length, bin_size):
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.spines['bottom'].set_position('zero')
    # configure spines of the graph; 'remove' top and right spines
    
    axs.spines['left'].set_bounds(0, max_y)
    axs.spines['bottom'].set_bounds(min_x, max_x)  # bottom line starts and end with length of genome
    axs.spines['left'].set_position(('data', min_x))

    # axis labeling
    #plt.xlabel("E. coli genomic coordinate (Mbp)")
    #plt.ylabel("Read count")
    axs.yaxis.set_label_coords(-0.02, 0.5)

    # set up xticks
    increments = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000]
    i = 0
    increment = increments[0]
    while (x_boundaries[1]-x_boundaries[0]) // increment > 10:
        i += 1
        increment = increments[i]
    
    x_axis_size = int(genome_length / bin_size)

    x_ticks = [r for r in range(0, x_axis_size, int(increment / bin_size)) if r >= (x_boundaries[0])/bin_size and r<=(x_boundaries[1])/bin_size]
    axs.set_xticks(x_ticks)
    x_tick_labels = [round(r, 2) for r in np.arange(x_boundaries[0]/1000, (x_boundaries[1]) / 1000, increment / 1000)]
    axs.set_xticklabels(x_tick_labels)

    # set limits of graph area - different from spine limits
    # leave extra space on top for info text, and space on bottom for the spacer marker
    plt.gca().set_xlim(left=-0.15*(max_x-min_x) + min_x, right=max_x)
    plt.gca().set_ylim(bottom=-0.12*max_y, top=1.20*max_y)

    return axs

def get_bins(filename, genome_length, bin_size):
    csv = np.genfromtxt(filename, delimiter=",", skip_header=1)
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

def plot_binned(filepath, run_information, yAxis_type, isPlasmid=False):  
    # Main graphing function for histogram binning
    # The yAxis_type can one of either:
    # 1. "raw" - the raw counts of reads in each bin
    # 2. "normalized" - by percentage of total reads
    # 3. "zoomed" - a normalized graph that is zoomed in at the y-axis to show low-frequence bins
    spacer_locations = [int(loc)+49 for loc in run_information['End of protospacer'].split()]
    genome_path = run_information['Target fasta file']
    if isPlasmid:
        genome_path = run_information['Second target fasta file']
    genome_length = len(SeqIO.read(Path(genome_path), 'fasta'))
    bin_scale = int(100000/genome_bin_size)
    bin_size = genome_bin_size
    if isPlasmid:
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
    max_x = x_boundaries[1] // genome_bin_size
    min_x = x_boundaries[0] // genome_bin_size

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
    setup_axes(axs, min_x, max_x, max_y, genome_length, bin_size)
    # 2 ticks on the y axis, one at 0 and on at max value of y (max reads)
    axs.set_yticks([0, max_y])
    axs.set_yticklabels([0, max_y_label])
    
    # 1. scatter plot points for protospacer markers
    for spacer_bin in spacer_bins:
        axs.scatter(spacer_bin, -0.03 * max_y, marker="^", c="#A32E79", s=33)
    # 2. main bar graph, width slightly >1 to avoid gaps between bars
    start_bin = x_boundaries[0]//bin_size
    end_bin = (x_boundaries[1]//bin_size)+1
    axs.bar(b[start_bin:end_bin], a2[start_bin:end_bin], color='#0D62AC', width=1.05)  

    # size of output figures
    fig.set_size_inches(fig_size_inches[0], fig_size_inches[1])

    # The text on the graph
    text_x = 5 * bin_scale
    text_y = 1.17 * max_y
    if yAxis_type == 'zoomed':
        text_x = 33 * bin_scale
        text_y = 0.8 * max_y

    # save figure
    sample = run_information['Sample']
    plt.savefig(os.path.join('outputs', 'plots', '{}_custom_bins_hist_{}.{}'.format(sample, yAxis_type, plots_filetype)), dpi=plots_dpi)
    plt.close()  # closes the matplotlib preview popup window


# PUT THIS SOMEWHERE IN THE TOP LEVEL OF THE REPO - eg. same level as requirements.txt, main.py, etc

#change these things - the meta is probably same as whatever the input file has
meta = {
    'Sample': 'A4574',
    'Genome': './bl21_de3_wt.fasta',
    'End of protospacer': '335151'
}
# path to the output csv of read locations
genome_reads_csv_path = '../genome_read_locations.csv'
# x boundaries for the graph, as from the genome refseq numbers
x_boundaries = [335000, 336000]
# raw, normalized, or zoomed
plot_type = 'zoomed'

# set the genome_bin_size and low_reads_cap_percent (for zoomed plots) in the parameters.py

# outputs will go in the normal output/plots folder with 'custom_bins_hist' in the name instead of 'genome'

#leave this alone
plot_binned(genome_reads_csv_path, meta, plot_type, False)