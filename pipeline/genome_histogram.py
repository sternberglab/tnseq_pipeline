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
exp_date = '191220'
percent = 0.50

# Run parameters for transposition distance histogram
query_length = 500
on_target_window = 100
plot_overlap = True

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
    genome_length = int(run_information['Genome Length'])
    desc = run_information['Information for graphs']
    psl = run_information['pCascade #']
    spacer_location = int(run_information['End of protospacer'])

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

def make_trans_dist_plot(fastaFile, run_information):
    # map spacer to refseq and determine query window
    

    code = run_information['Sample']
    qScore = run_information['Qscore Threshold']
    genome_length = int(run_information['Genome Length'])
    desc = run_information['Information for graphs']
    psl = run_information['pCascade #']
    spacer_location = int(run_information['End of protospacer'])
    description = run_information['Description']
    spacer = run_information['Spacer']
    
    genome = SeqIO.read(Path(genome_path), "fasta")
    refseq = genome.seq.upper()
    if run_information['Target direction'].lower() == 'rv':
        refseq = refseq.reverse_complement()

    if refseq.find(spacer) >= 0:
        spacer_end = refseq.find(spacer) + 32
        query = refseq[spacer_end-90:spacer_end+query_length]  # '-90' accounts for the -50bp of the on-target later
        query_rc = query.reverse_complement()
        spacer_end = 90  # resets spacer end index to middle of query window (no longer using full refseq)
    else:
        print("ERROR - Spacer not found within RefSeq")
        return
    total = 0  # counts total reads in fastq file
    out_list_all = []  # list holding common trans_dist
    out_list_rl = []  # list holding indv RL trans_dist values
    out_list_lr = []  # list holding indv LR trans_dist values
    # these lists are longer than query_length to hold negative values of trans_dist
    # the output excel cuts of list using query_length so those values will not show
    out_tally_rl = [0] * (query_length+spacer_end+20)  # list tallying freq of tran_dist for RL
    out_tally_lr = [0] * (query_length+spacer_end+20)  # list tallying freq of tran_dist for LR
    example_reads_rl = ['X'] * (query_length+spacer_end+20)  # to hold example reads mapping to each trans_dist for RL
    example_reads_lr = ['X'] * (query_length+spacer_end+20)  # to hold example reads mapping to each trans_dist for LR
    for record in SeqIO.parse(fastaFile, Path(fastaFile).suffix[1:]):  # loop through all reads in fastq file
        total += 1
        new_seq = record.seq  # artifact from v1 of code
        if query.find(new_seq) >= 0:  # corresponds to RL
            trans_dist = query.find(new_seq) + 17 - spacer_end  # distance in bp from end of protospacer
            out_list_all.append(trans_dist)  # append to holding lists for processing later
            out_list_rl.append(trans_dist)
            if out_tally_rl[trans_dist] == 0:  # add read to example list if this is the first occurrence
                example_reads_rl[trans_dist] = new_seq
            out_tally_rl[trans_dist] += 1  # count into tally list
        elif query_rc.find(new_seq) >= 0:  # corresponds to LR
            trans_dist = len(query) - query_rc.find(new_seq) - 17 + 5 - spacer_end  # dist in bp from end of protospacer
            out_list_all.append(trans_dist)  # append to tally lists for processing later
            out_list_lr.append(trans_dist)
            if out_tally_lr[trans_dist] == 0:  # add read to example list if this is the first occurrence
                example_reads_lr[trans_dist] = new_seq
            out_tally_lr[trans_dist] += 1  # count into tally list
    #  determine most frequent trans_dist
    out_tally_all = [0] * query_length
    for i in range(0, len(out_tally_all)):
        out_tally_all[i] = out_tally_rl[i] + out_tally_lr[i]
    for x, y in enumerate(out_tally_all):
        if y == max(out_tally_all):
            main_site = x + spacer_end  # remember to convert dist to site of integration

    # define on target window
    on_target_lower = main_site - int(on_target_window/2)
    on_target_upper = main_site + int(on_target_window/2)
    # move any trans_dist within this window into a final holding list and clears old holding list
    final_list_rl = []
    for dist in out_list_rl:
        if on_target_lower <= (dist + spacer_end) <= on_target_upper:  # convert dist to site of integration
            final_list_rl.append(dist)

    final_list_lr = []
    for dist in out_list_lr:
        if on_target_lower <= (dist + spacer_end) <= on_target_upper:  # convert dist to site of integration
            final_list_lr.append(dist)

    # determine on target frequency
    on_target_total = len(final_list_rl) + len(final_list_lr)
    off_target = total - on_target_total

    # determine top 3 most common trans_dist for highlight box
    # for combined RL and LR
    indices = []  # for zipping with out_tally lists
    for i in range(0, query_length):
        indices.append(i)
    top_3 = heapq.nlargest(3, zip(out_tally_all, indices))  # exists as a list of smaller 2-item lists

    x_axis = []  # artificial x-axis
    for i in range(20, 61):
        x_axis.append(i)
    y_rl = out_tally_rl[20:61]
    y_lr = out_tally_lr[20:61]

    max_y = max(max(y_rl), max(y_lr))
    fig, axs = plt.subplots(1, 1, tight_layout=True)
    title = fig.suptitle("{} - {} / On-target = {}% / Bias = {} :1".format(
        code, description, round(100*on_target_total/total, 1),
                         round(len(final_list_rl)/(len(final_list_lr)+0.00000001), 2)))
    title.set_y(0.9)
    # LR graph is colorless with a border
    axs.bar(x_axis, y_lr, color='none', edgecolor='#153C6B', linewidth=1.0, width=1.01, zorder=1)
    # RL graph is blue with no border (behind bordered RL)
    axs.bar(x_axis, y_rl, color='#83B0DD', edgecolor='#83B0DD', linewidth=1.0, width=1.01, zorder=0)

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
    # plt.xlabel("Distance from target site (bp)")
    # plt.ylabel("Read count")
    # plt.gca().set_xlim(left=35, right=60)
    # plt.gca().set_ylim(bottom=0, top=(1.3*max_y))

    # plt.savefig('test.svg', dpi=250)
    run_prefix = "{}_Q{}".format(code, qScore)
    plt.savefig(output_path(os.path.join('plots', '{}_trans_dist_hist.svg'.format(run_prefix))), dpi=500)
    plt.close()


def make_plots(csvFile, fastaFile):
    start = time.perf_counter()
    print("Making genome mapping histograms and transposition distance plot...")
    # get code and q score from the filename, ex. "A4632_Q20_unique_reads.csv"
    code, qScore = Path(csvFile).stem.split('_')[:2]
    qScore = qScore[1:]
    code_info = None
    code_logs = None
    # read in the info csv files
    with open(Path(info_file), 'r', encoding='utf-8-sig') as info_csv:
        reader = csv.DictReader(info_csv)
        for row in reader:
            if row['Sample'] == code:
                code_info = row
    code_logs = get_log_entry(code, qScore)
    run_information = {**code_logs, **code_info}
    print("Got the meta information about this run, creating the genome-mapping histograms...")
    plot_binned(csvFile, run_information, 'raw')
    plot_binned(csvFile, run_information, 'normalized')
    plot_binned(csvFile, run_information, 'zoomed')
    print("Doing the transposition distance mapping...")
    make_trans_dist_plot(fastaFile, run_information)

    elapsed_time = time.perf_counter() - start
    print("Finished all plotting in {}".format(elapsed_time))
