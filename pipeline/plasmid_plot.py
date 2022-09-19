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

from parameters import info_file, plots_filetype, plots_dpi, fig_size_inches
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

plasmid_bin_size = 20
xtick_increment = 1000

def get_bins(filename, genome_length):
    csv = np.genfromtxt(filename, delimiter=",", skip_header=1)
    number_of_bins = (genome_length // plasmid_bin_size) + 1
    bin_numbers = np.array(range(1, number_of_bins + 1))
    bin_values = []
    for i in bin_numbers:
        min_val = (i - 1) * plasmid_bin_size
        max_val = i * plasmid_bin_size
        counts_in_bin = [entry[1] for entry in csv if min_val <= entry[0] < max_val]
        total_for_bin = sum(counts_in_bin)
        bin_values.append(total_for_bin)
    np_bin_values = np.asarray(bin_values)
    return bin_numbers, np_bin_values

def setup_axes(axs, max_x, max_y, genome_length):
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
    
    
    x_axis_size = int(genome_length / plasmid_bin_size)+1

    axs.set_xticks(range(0, x_axis_size, int(xtick_increment / plasmid_bin_size)))
    axs.set_xticklabels(range(0, int(genome_length / xtick_increment)+1))

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
    
    genome_path = run_information['Second target fasta file']
    genome_length = len(SeqIO.read(Path(genome_path), 'fasta'))
        
    bin_scale = int(genome_length // 500)
    
    # determine which bin the spacer lies in
    spacer_bins = []
    if not isPlasmid:
        spacer_bins = [int(spacer_location / plasmid_bin_size) + 0.5 for spacer_location in spacer_locations]

    # get bins and counts for the histogram
    try:
        b, a2 = get_bins(filepath, genome_length)
    except FileNotFoundError as e:
        return

    total_reads = int(sum(a2))

    max_x = len(a2)

    # For normalized and zoomed graphs, normalized based on total reads as a percentage
    max_y = max_y_label = None
    if yAxis_type == 'normalized':
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
    setup_axes(axs, max_x, max_y, genome_length)
    
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
    graph_name = "second_target"
    plot_path = output_path(os.path.join('plots', f'{sample}_{graph_name}_hist_{yAxis_type}.{plots_filetype}'))
    plt.savefig(plot_path, dpi=plots_dpi)
    plt.close()  # closes the matplotlib preview popup window

def plot_section(csvFile, meta_info, spacerSeq, name):
    genome_path = meta_info['Second target fasta file']
    genome = SeqIO.read(Path(genome_path), 'fasta')

    spacerEnd = genome.seq.upper().find(spacerSeq) + len(spacerSeq)

    reads = []
    try:
        with open(csvFile, encoding='utf-8-sig') as openFile:
            reader = csv.DictReader(openFile)
            reads = list(reader)
    except FileNotFoundError as e:
        return
    
    x_axis = list(range(40,61))
    y_vals = []
    for i in range(spacerEnd + 40, spacerEnd + 61):
        match = next((item for item in reads if int(item['position']) == i), None)
        if match:
            y_vals.append(int(match['reads']))
        else:
            y_vals.append(0)

    max_y = max(y_vals)

    fig, axs = plt.subplots(1, 1, tight_layout=True)
    # Blue with no border for total numbers
    axs.bar(x_axis, y_vals, color='#83B0DD', edgecolor='#83B0DD', linewidth=1.0, width=1.01, zorder=0)

    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    # axs.spines['left'].set_visible(False)
    axs.spines['bottom'].set_position('zero')
    axs.spines['left'].set_bounds(0, max_y)
    axs.set_xticks([40, 42, 45, 50, 55, 60])
    axs.set_xticklabels([40, 0, 45, 50, 55, 60])
    axs.set_yticks([0, max_y])
    axs.set_yticklabels([0, max_y])
    axs.set_xlim(left=42, right=58) ## Change window here
    axs.set_ylim(bottom=0, top=1.25 * (max_y))
    axs.set(xlabel="Distance from target site (bp)", ylabel="Read count")
    axs.yaxis.set_label_coords(-0.05, 0.4)

    fig.set_size_inches(5, 4.2)
    plot_path = output_path(os.path.join('plots', f'{meta_info["Sample"]}_second_target_dist_{name}.{plots_filetype}'))
    plt.savefig(plot_path, dpi=plots_dpi)
    plt.close()


def plot_plasmid(csvFile, meta_info):
    CRISPR_array_spacer = meta_info['CRISPR Array Sequence']
    donor_spacer = meta_info['Donor sequence']
    start = time.perf_counter()
    print("Making plasmid mapping histograms...")

    print("Got the meta information about this run, creating the plasmid-mapping histograms...")
    plot_binned(csvFile, meta_info, 'raw')
    plot_binned(csvFile, meta_info, 'normalized')
    # plot sections for the crispr array region and donor region
    if CRISPR_array_spacer:
        plot_section(csvFile, meta_info, CRISPR_array_spacer, "CRISPR_seq")
    if donor_spacer:
        plot_section(csvFile, meta_info, donor_spacer, "donor")
    elapsed_time = round(time.perf_counter() - start, 2)
    print("Finished plasmid mapping plotting in {} seconds".format(elapsed_time))
    return meta_info
