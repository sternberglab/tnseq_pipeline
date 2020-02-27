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

from parameters import info_file, genome_path
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

# Run parameters for genome wide histograms
bin_size = 5000  # in bp
bin_scale = int(100000/bin_size)
percent = 0.50

def get_bins(filename, bin_size, genome_length):
    csv = np.genfromtxt(filename, delimiter=",", skip_header=1)
    number_of_bins = (genome_length // bin_size) + 1
    bin_numbers = np.array(range(1, number_of_bins + 1))
    bin_values = []
    for i in bin_numbers:
        min_val = (i - 1) * bin_size
        max_val = i * bin_size
        counts_in_bin = [entry[1] for entry in csv if min_val <= entry[0] < max_val]
        total_for_bin = sum(counts_in_bin)
        bin_values.append(total_for_bin)
    np_bin_values = np.asarray(bin_values)
    return bin_numbers, np_bin_values

def setup_axes(axs, max_x, max_y):
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
    axs.set_xticks(range(0, 50*bin_scale, 5*bin_scale))

    axs.set_xticklabels(np.arange(0, 5, 0.5))

    # set limits of graph area - different from spine limits
    # leave extra space on top for info text, and space on bottom for the spacer marker
    plt.gca().set_xlim(left=-0.15*max_x, right=max_x)
    plt.gca().set_ylim(bottom=-0.12*max_y, top=1.20*max_y)

    return axs

def plot_binned(filepath, run_information, yAxis_type):  
    # Main graphing function for histogram binning
    # The yAxis_type can one of either:
    # 1. "raw" - the raw counts of reads in each bin
    # 2. "normalized" - by percentage of total reads
    # 3. "zoomed" - a normalized graph that is zoomed in at the y-axis to show low-frequence bins
    code = run_information['Sample']
    qScore = run_information['Qscore Threshold']
    desc = run_information['Information for graphs']
    psl = run_information['pCascade #']
    exp_date = run_information['Experiment date']
    spacer_location = int(run_information['End of protospacer'])

    genome_length = len(SeqIO.read(Path(genome_path), 'fasta'))

    # determine which bin the spacer lies in
    spacer_bin = int(spacer_location/bin_size) + 0.5  

    # get bins and counts for the histogram
    b, a2 = get_bins(filepath, bin_size, genome_length)

    total_reads = int(sum(a2))

    max_x = len(a2)

    # For normalized and zoomed graphs, normalized based on total reads as a percentage
    max_y = max_y_label = None
    if yAxis_type == 'zoomed':
        max_y = total_reads * percent * 0.01
        max_y_label = percent
    elif yAxis_type == 'normalized':
        a2 = (100*a2)/total_reads
        max_y = 100
        max_y_label = '100%'
    elif yAxis_type == 'raw':
        max_y = int(max(a2))
        max_y_label = str(max_y)
    else:
        print('Must give a yAxis type or raw, normalized, or zoomed')

    # set up figure
    fig, axs = plt.subplots(1, 1, tight_layout=True)
    setup_axes(axs, max_x, max_y)
    
    # 2 ticks on the y axis, one at 0 and on at max value of y (max reads)
    axs.set_yticks([0, max_y])
    axs.set_yticklabels([0, max_y_label])
    
    # 1 scatter plot point for protospacer marker
    axs.scatter(spacer_bin, -0.03 * max_y, marker="^", c="#A32E79", s=33)
    # 2. main bar graph, width slightly >1 to avoid gaps between bars
    axs.bar(b, a2, color='#0D62AC', width=1.05)  

    # size of output figures
    fig.set_size_inches(8.9, 3.1)

    # The text on the graph
    text_x = 5 * bin_scale
    text_y = 1.17 * max_y
    if yAxis_type == 'zoomed':
        text_x = 33 * bin_scale
        text_y = 0.8 * max_y
    plt.text(text_x, text_y, "{}-{}\n{}\ngRNA Plasmid = {}\nTotal Reads = {}"
        .format(code, exp_date, desc, psl, total_reads))

    # save figure
    run_prefix = "{}_Q{}".format(code, qScore)
    plt.savefig(output_path(os.path.join('plots', '{}_genome_hist_{}.svg'.format(run_prefix, yAxis_type))), dpi=500)
    plt.close()  # closes the matplotlib preview popup window


def make_genome_plots(csvFile, meta_info):
    start = time.perf_counter()
    print("Making genome mapping histograms...")
    # get code and q score from the filename, ex. "A4632_Q20_unique_reads.csv"
    code, qScore = Path(csvFile).stem.split('_')[:2]
    qScore = qScore[1:]

    code_info = meta_info
    code_logs = get_log_entry(code, qScore)
    run_information = {**code_logs, **code_info}
    print("Got the meta information about this run, creating the genome-mapping histograms...")
    plot_binned(csvFile, run_information, 'raw')
    plot_binned(csvFile, run_information, 'normalized')
    plot_binned(csvFile, run_information, 'zoomed')
    elapsed_time = round(time.perf_counter() - start, 2)
    print("Finished genome mapping plotting in {} seconds".format(elapsed_time))
    return run_information
