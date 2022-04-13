import os
import fnmatch
from Bio import SeqIO
import xlsxwriter
from pathlib import Path
import datetime as dt
import csv
import matplotlib.pyplot as plt
import time

from parameters import on_target_window, plots_filetype, plots_dpi
from .utils import output_path

from types import SimpleNamespace

# The plot boundaries, as distance from end of protospacer
PLOT_BOUNDARIES = [10, 61]


def make_trans_dist_plot(readsCsv, run_information):
    print("Doing the transposition distance mapping...")
    start = time.perf_counter()

    # map spacer to refseq and determine query window
    code = run_information.get('Sample')
    genome_path = run_information.get('Genome')
    genome_length = len(SeqIO.read(Path(genome_path), 'fasta'))
    desc = run_information.get('Information for graphs', None)
    psl = run_information.get('pCascade #', None)
    description = run_information.get('Description', None)
    spacer = run_information.get('Spacer').upper()
    exp_date = run_information.get('Experiment date', None)

    genome = SeqIO.read(Path(genome_path), "fasta")
    refseq = genome.seq.upper()
    is_reverse = False
    if refseq.find(spacer) >= 0:
        spacer_end = refseq.find(spacer) + len(spacer)
        target_site = spacer_end + 50
    elif refseq.reverse_complement().find(spacer) >= 0:
        spacer_end = len(refseq) - (refseq.reverse_complement().find(spacer) + len(spacer))
        target_site = spacer_end - 50
        is_reverse = True
    else:
        print("ERROR - Spacer not found within RefSeq")
        return

    all_reads = []
    with open(readsCsv, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        all_reads = list([row for row in reader])

    x_axis = []  # artificial x-axis
    for i in range(PLOT_BOUNDARIES[0], PLOT_BOUNDARIES[1]+1):
        x_axis.append(i)
    y_rl = []
    y_lr = []
    total = sum([int(row["reads"]) for row in all_reads])
    reads_in_window = [row for row in all_reads if abs(int(row["position"]) - target_site) < (int(on_target_window) / 2)]
    for i in range(PLOT_BOUNDARIES[0], PLOT_BOUNDARIES[1]+1):
        if is_reverse:
            reads = next(iter([row for row in all_reads if int(row["position"]) == spacer_end - i]), None)
        else:
            reads = next(iter([row for row in all_reads if int(row["position"]) == spacer_end + i]), None)
        if reads:
            y_rl.append(int(float(reads["RL"])))
            y_lr.append(int(float(reads["LR"])))
        else:
            y_rl.append(0)
            y_lr.append(0)
    on_target_total = sum([int(row["reads"]) for row in reads_in_window])
    off_target_reads = total - on_target_total

    max_y = max(max(y_rl), max(y_lr))

    tick_distance = int((PLOT_BOUNDARIES[1] - PLOT_BOUNDARIES[0])/5)
    x_ticks = range(PLOT_BOUNDARIES[0], PLOT_BOUNDARIES[1] + 2, tick_distance)
    for overlap in [True, False]:
        if overlap:
            fig, axs = plt.subplots(1, 1, tight_layout=True)
            title = fig.suptitle("{} - {} / On-target = {}% / Bias = {} :1".format(
                code, description, round(100*on_target_total/total, 1),
                                 round(sum(y_rl)/(sum(y_lr)+0.00000001), 2)))
            title.set_y(0.9)

            # LR graph is colorless with a border
            axs.bar(x_axis, y_lr, color='none', edgecolor='#153C6B', linewidth=1.0, width=1.01, zorder=1)
            # RL graph is blue with no border (behind bordered RL)
            axs.bar(x_axis, y_rl, color='#83B0DD', edgecolor='#83B0DD', linewidth=1.0, width=1.01, zorder=0)

            axs.spines['top'].set_visible(False)
            axs.spines['right'].set_visible(False)
            axs.spines['bottom'].set_position('zero')
            axs.spines['left'].set_bounds(0, max_y)
            axs.set_xticks(x_ticks)
            axs.set_xticklabels(x_ticks)
            axs.set_yticks([0, max_y])
            axs.set_yticklabels([0, max_y])
            axs.set_xlim(left=PLOT_BOUNDARIES[0], right=PLOT_BOUNDARIES[1])
            axs.set_ylim(bottom=0, top=1.25 * (max_y))
            axs.set(xlabel="Distance from target site (bp)", ylabel="Read count")
            axs.yaxis.set_label_coords(-0.05, 0.4)

            fig.set_size_inches(5, 4.2)
            plot_name = "{}_trans_dist_hist_overlap_v2.svg".format(code)
        else:
            fig, axs = plt.subplots(1, 2)
            fig.tight_layout(rect=[0.15, 0.1, 1, 0.9])
            title = fig.suptitle("{} - {}\nOn-target frequency: {}%\nOrientation bias (R->L:L->R): {}:1"
                                 .format(code, description, round(100*on_target_total/total, 1),
                                 round(sum(y_rl)/(sum(y_lr)+0.00000001), 2)))
            title.set_y(0.88)
            # first graph on the left in red
            axs[0].bar(x_axis, y_rl, color='tab:orange', width=1.0)
            axs[0].set_title("R->L Integration Events")
            # second graph on the right in blue
            axs[1].bar(x_axis, y_lr, color='tab:blue', width=1.0)
            axs[1].set_title("L->R Integration Events")
            fig.subplots_adjust(wspace=0.7)

            for axs in axs.flat:
                axs.spines['top'].set_visible(False)
                axs.spines['right'].set_visible(False)
                axs.spines['bottom'].set_position('zero')
                axs.spines['left'].set_bounds(0, max_y)
                axs.set_xticks(x_ticks)
                axs.set_xticklabels(x_ticks)
                axs.set_yticks([0, max_y])
                axs.set_yticklabels([0, max_y])
                axs.set_xlim(left=PLOT_BOUNDARIES[0], right=PLOT_BOUNDARIES[1])
                axs.set_ylim(bottom=0, top=1.05*(max_y))
                axs.set(xlabel="Distance from target site (bp)", ylabel="Read count")
                axs.yaxis.set_label_coords(-0.1,0.5)

            fig.set_size_inches(6, 4.2)
            fig.subplots_adjust(top=0.65)
            plot_name = "{}_trans_dist_hist_no_overlap_v2.{}".format(code, plots_filetype)
        # save the plot and close it
        plt.savefig(output_path(os.path.join('plots', plot_name)), dpi=plots_dpi)
        plt.close()

    # Return information to append to the overall output logs for this run (Sample code plus Qscore combo)
    # create_output_excel(SimpleNamespace(**locals()))
    elapsed_time = round(time.perf_counter() - start, 2)
    print("Finished transposition distance plotting and output in {} seconds".format(elapsed_time))
    return