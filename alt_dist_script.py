import os
import fnmatch
from Bio import SeqIO
import xlsxwriter
from pathlib import Path
import datetime as dt
import csv
import heapq
import matplotlib.pyplot as plt
import time

from parameters import query_length, on_target_window, plots_filetype, plots_dpi

from types import SimpleNamespace

# Make an output CSV with Position, RL, LR, and Total columns for this to read in
# Need to manually copy over on-target % since the output values for the CSV are 
# just for the query window (500bp)
csvFile = '/Users/Chris_Acree/Documents/lab/working/927-913_combined.csv'
on_target1 = 13137
off_target1 = 15125
on_target2 = 9854
off_target2 = 12390

# To adjust the plotting window, change the 'axs.set_xlim' parameters for right and left
xlim_left = 28
xlim_right = 52

def make_trans_dist_plot():
    run_prefix = Path(csvFile).stem.split('_')[0]

    with open(csvFile, 'r', encoding='utf-8') as opened:
        reader = csv.DictReader(opened)
        rows = list(reader)
        y_total = [int(row['Total']) for row in rows]
        out_tally_rl = [int(row['RL']) for row in rows]
        out_tally_lr = [int(row['LR']) for row in rows]


    x_axis = []  # artificial x-axis
    for i in range(0, 61):
        x_axis.append(i)
    y_rl = out_tally_rl[0:61]
    y_lr = out_tally_lr[0:61]

    max_y = max(max(y_rl), max(y_lr))
    
    on_target_perc = round(100*(on_target1+on_target2)/(on_target1+on_target2+off_target1+off_target2), 2)
    bias = round(sum(out_tally_rl)/sum(out_tally_lr), 2)

    xticks = range(0,61,5)

    for overlap in [True, False]:
        if overlap:
            fig, axs = plt.subplots(1, 1, tight_layout=True)
            title = fig.suptitle("{} / On-target = {}% / Bias = {} : 1".format(
                run_prefix, on_target_perc, bias))

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
            axs.set_xticks(xticks)
            axs.set_xticklabels(xticks)
            axs.set_yticks([0, max_y])
            axs.set_yticklabels([0, max_y])
            axs.set_xlim(left=xlim_left, right=xlim_right) ## Change window here
            axs.set_ylim(bottom=0, top=1.25 * (max_y))
            axs.set(xlabel="Distance from target site (bp)", ylabel="Read count")
            axs.yaxis.set_label_coords(-0.05, 0.4)

            fig.set_size_inches(5, 4.2)
            plot_name = "{}_trans_dist_hist_overlap.{}".format(run_prefix, plots_filetype)
        else:
            fig, axs = plt.subplots(1, 2)
            fig. tight_layout(rect=[0.15, 0.1, 1, 0.9])
            title = fig.suptitle("{} \nOn-target frequency: {}%\nOrientation bias (R->L:L->R): {}:1"
                                 .format(run_prefix, on_target_perc, bias))
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
                # ax.spines['bottom'].set_visible(False)
                # axs.spines['left'].set_visible(False)
                axs.spines['bottom'].set_position('zero')
                axs.spines['left'].set_bounds(0, max_y)
                axs.set_xticks(xticks)
                axs.set_xticklabels(xticks)
                axs.set_yticks([0, max_y])
                axs.set_yticklabels([0, max_y])
                axs.set_xlim(left=xlim_left, right=xlim_right) ## Change window here
                axs.set_ylim(bottom=0, top=1.05*(max_y))
                axs.set(xlabel="Distance from target site (bp)", ylabel="Read count")
                axs.yaxis.set_label_coords(-0.1,0.5)

            fig.set_size_inches(6, 4.2)
            fig.subplots_adjust(top=0.65)
            plot_name = "{}_trans_dist_hist_no_overlap.{}".format(run_prefix, plots_filetype)
        # save the plot and close it
        plt.savefig(plot_name, dpi=plots_dpi)
        plt.close()

    # Return information to append to the overall output logs for this run (Sample code plus Qscore combo)
    print("Finished transposition distance plotting and output")
    return

make_trans_dist_plot()