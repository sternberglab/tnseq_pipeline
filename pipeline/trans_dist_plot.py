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

'''
def create_output_excel(v):
    excel_path = Path(output_path('{}_trans_dist_mmei_output.xlsx'.format(v.run_prefix)))

    log = xlsxwriter.Workbook(excel_path)
    bold = log.add_format({'bold': True})
    upsizebold = log.add_format()
    upsizebold.set_font_size(16)
    upsizebold.set_bold()
    percentage_format = log.add_format()
    percentage_format.set_num_format('0.00%')
    deci3_format = log.add_format()
    deci3_format.set_num_format('0.000')
    red_font = log.add_format({'font_color': 'red'})
    blue_font = log.add_format({'font_color': 'blue'})
    red_percent = log.add_format()
    red_percent.set_num_format('0.00%')
    red_percent.set_font_color('red')
    blue_percent = log.add_format()
    blue_percent.set_num_format('0.00%')
    blue_percent.set_font_color('blue')
    red_deci3 = log.add_format()
    red_deci3.set_num_format('0.000')
    red_deci3.set_font_color('red')
    blue_deci3 = log.add_format()
    blue_deci3.set_num_format('0.000')
    blue_deci3.set_font_color('blue')

    # make landing info sheet containing all information shared by sheets in the excel file
    infosheet = log.add_worksheet("General Info")
    infosheet.write(0, 0, "NGS MmeI Analysis of Transposition Distance ", upsizebold)
    infosheet.set_column(0, 0, 35)
    infosheet.write(2, 0, "NextSeq/NGS Data Collection Date", bold)
    infosheet.write(2, 1, v.exp_date)

    infosheet.write(3, 0, "Python Data Analysis Date", bold)
    infosheet.write(3, 1, str(dt.datetime.now())[0:10])

    infosheet.write(4, 0, "Username/Initials", bold)
    infosheet.write(4, 1, "")

    infosheet.write(5, 0, "Python Code Used", bold)
    infosheet.write(5, 1, "")

    infosheet.write(6, 0, "Notes", bold)
    infosheet.write(7, 0, '"Query Window" - larger window for query distances (not the same as on-target window'
                          ', same strand as the BL21 RefSeq')
    infosheet.write(8, 0, '"Genomic Base" - Base directly 5\' of integration site, on the SAME strand as protospacer')
    infosheet.write(9, 0, '"Example Reads" - Reads here are trimmed to 17bp and the reverse complement of reads in'
                          'the Geneious fastq output', )

    logsheet = log.add_worksheet(v.code)
    logsheet.set_column(0, 0, 24)
    logsheet.set_column(1, 1, 20)
    logsheet.set_column(2, 6, 17)
    logsheet.set_column(7, 12, 19)
    logsheet.write(3, 3, " ")  # for clearer aesthetic

    logsheet.write(0, 0, "Sample ID", bold)
    logsheet.write(0, 1, v.code)

    logsheet.write(1, 0, "Description", bold)
    logsheet.write(1, 1, v.description)

    logsheet.write(2, 0, "Target Location", bold)
    if v.direction == 'fw':
        logsheet.write(2, 1, "5' of Integration Site")
    else:
        logsheet.write(2, 1, "3' of Integration Site (RevCom)")

    logsheet.write(3, 0, "Query Window", bold)
    if v.direction == 'fw':
        logsheet.write(3, 1, str(v.query))
    if v.direction == 'rv':
        logsheet.write(3, 1, str(v.query_rc))

    logsheet.write(4, 0, "Plasmid encoding gRNA", bold)
    logsheet.write(4, 1, v.psl)

    logsheet.write(5, 0, "Protospacer", bold)
    logsheet.write(5, 1, v.spacer)

    logsheet.write(6, 0, "Total Reads", bold)
    logsheet.write(6, 1, v.total)

    logsheet.write(7, 0, "On Target Reads", bold)
    logsheet.write(7, 1, v.on_target_total)
    logsheet.write(7, 2, v.on_target_total / v.total, percentage_format)

    logsheet.write(8, 0, "Off Target Reads", bold)
    logsheet.write(8, 1, v.off_target)
    logsheet.write(8, 2, v.off_target / v.total, percentage_format)

    logsheet.write(9, 0, "On Target Reads in RL Orientation", bold)
    logsheet.write(9, 1, len(v.final_list_rl))
    logsheet.write(9, 2, len(v.final_list_rl) / v.total, percentage_format)

    logsheet.write(10, 0, "On Target Reads in LR Orientation", bold)
    logsheet.write(10, 1, len(v.final_list_lr))
    logsheet.write(10, 2, len(v.final_list_lr) / v.total, percentage_format)

    logsheet.write(11, 1, "Protospacer-Transposon Distance", bold)
    for i in range(0, query_length):
        logsheet.write(i + 12, 1, i)

    logsheet.write(11, 0, "Genomic Base", bold)
    for i in range(-1, query_length - 1):
        logsheet.write(i + 13, 0, v.query[i+v.spacer_end])  # shift back 1 to get the base right before transposition

    logsheet.write(11, 2, "Number of Reads (RL)", bold)
    for i in range(0, query_length):
        logsheet.write(i + 12, 2, v.out_tally_rl[i], red_font)
    logsheet.write(11, 3, "% of Total Reads (RL)", bold)
    for i in range(0, query_length):
        logsheet.write(i + 12, 3, v.out_tally_rl[i] / v.total, red_percent)
    logsheet.write(11, 4, "Normalized Read Count (RL)", bold)
    if max(v.out_tally_rl) > 0:
        for i in range(0, query_length):
            logsheet.write(i + 12, 4, v.out_tally_rl[i] / max(v.out_tally_rl), red_deci3)

    logsheet.write(11, 5, "Number of Reads (LR)", bold)
    for i in range(0, query_length):
        logsheet.write(i + 12, 5, v.out_tally_lr[i], blue_font)
    logsheet.write(11, 6, "% of Total Reads (LR)", bold)
    for i in range(0, query_length):
        logsheet.write(i + 12, 6, v.out_tally_lr[i] / v.total, blue_percent)
    logsheet.write(11, 7, "Normalized Read Count (LR)", bold)
    if max(v.out_tally_lr) > 0:
        for i in range(0, query_length):
            logsheet.write(i + 12, 7, v.out_tally_lr[i] / max(v.out_tally_lr), blue_deci3)

    logsheet.write(11, 8, "Number of Reads (Combined)", bold)
    for i in range(0, query_length):
        logsheet.write(i + 12, 8, v.out_tally_all[i])
    logsheet.write(11, 9, "% of Total Reads (Combined)", bold)
    for i in range(0, query_length):
        logsheet.write(i + 12, 9, v.out_tally_all[i] / v.total, percentage_format)
    logsheet.write(11, 10, "Normalized Read Count (Combined)", bold)
    if max(v.out_tally_all) > 0:
        for i in range(0, query_length):
            logsheet.write(i + 12, 10, v.out_tally_all[i] / max(v.out_tally_all), deci3_format)

    logsheet.write(11, 11, "Example Reads RL", bold)
    for i in range(0, query_length):
        logsheet.write(i + 12, 11, str(v.example_reads_rl[i]), red_font)

    logsheet.write(11, 12, "Example Reads LR", bold)
    for i in range(0, query_length):
        logsheet.write(i + 12, 12, str(v.example_reads_lr[i]), blue_font)

    # determine top 3 most common trans_dist for highlight box
    # for combined RL and LR
    indices = []  # for zipping with out_tally lists
    for i in range(0, query_length):
        indices.append(i)
    top_3 = heapq.nlargest(3, zip(v.out_tally_all, indices))  # exists as a list of smaller 2-item lists

    # 'highlight box', take from top_3 list determined above
    logsheet.write(0, 4, 'Most Frequent Transposition Distances (bp)', bold)

    logsheet.write(1, 4, top_3[0][1])
    logsheet.write(1, 5, top_3[0][0] / v.total, percentage_format)
    logsheet.write(1, 6, v.out_tally_rl[top_3[0][1]] / v.total, red_percent)
    logsheet.write(1, 7, v.out_tally_lr[top_3[0][1]] / v.total, blue_percent)

    logsheet.write(2, 4, top_3[1][1])
    logsheet.write(2, 5, top_3[1][0] / v.total, percentage_format)
    logsheet.write(2, 6, v.out_tally_rl[top_3[1][1]] / v.total, red_percent)
    logsheet.write(2, 7, v.out_tally_lr[top_3[1][1]] / v.total, blue_percent)

    logsheet.write(3, 4, top_3[2][1])
    logsheet.write(3, 5, top_3[2][0] / v.total, percentage_format)
    logsheet.write(3, 6, v.out_tally_rl[top_3[2][1]] / v.total, red_percent)
    logsheet.write(3, 7, v.out_tally_lr[top_3[2][1]] / v.total, blue_percent)

    logsheet.write(4, 4, 'On Target Frequency', bold)
    logsheet.write(4, 5, v.on_target_total / v.total, percentage_format)

    logsheet.write(5, 4, 'Orientation Bias (R->L:L->R)', bold)
    logsheet.write(5, 5, '{} : 1'.format(round(len(v.final_list_rl)/(len(v.final_list_lr)+0.00000001), 2)))  # in case LR is 0
    log.close()
    return
'''

def make_trans_dist_plot(readsCsv, run_information):
    print("Doing the transposition distance mapping...")
    start = time.perf_counter()

    # map spacer to refseq and determine query window
    code = run_information['Sample']
    genome_path = run_information['Genome']
    genome_length = len(SeqIO.read(Path(genome_path), 'fasta'))
    desc = run_information['Information for graphs']
    psl = run_information['pCascade #']
    description = run_information['Description']
    spacer = run_information['Spacer'].upper()
    exp_date = run_information['Experiment date']

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
    for i in range(20, 61):
        x_axis.append(i)
    y_rl = []
    y_lr = []
    total = sum([int(row["reads"]) for row in all_reads])
    reads_in_window = [row for row in all_reads if abs(int(row["position"]) - target_site) < (int(on_target_window) / 2)]
    for i in range(20, 61):
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
            plot_name = "{}_trans_dist_hist_overlap_v2.svg".format(code)
        else:
            fig, axs = plt.subplots(1, 2)
            fig. tight_layout(rect=[0.15, 0.1, 1, 0.9])
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
                # ax.spines['bottom'].set_visible(False)
                # axs.spines['left'].set_visible(False)
                axs.spines['bottom'].set_position('zero')
                axs.spines['left'].set_bounds(0, max_y)
                axs.set_xticks([40,45,50,55,60])
                axs.set_xticklabels([40,45,50,55,60])
                axs.set_yticks([0, max_y])
                axs.set_yticklabels([0, max_y])
                axs.set_xlim(left=39, right=61) ## Change window here
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