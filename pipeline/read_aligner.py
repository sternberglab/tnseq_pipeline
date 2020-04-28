import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import pandas as pd
import sys
import csv
import multiprocessing
import shutil
import numpy as np
import subprocess
import platform
import time

from .utils import inter_path, output_path
from parameters import fingerprint_length as map_length, transposon_site_duplication_length as TSD

def find_alignments(read_sequences_fasta, genome_fasta, output_filename):
    print("Find alignments for {}".format(output_filename))

    genome_fasta_path = Path(genome_fasta)
    if not genome_fasta_path.exists():
        print("Source sequence file not found at {}".format(genome_fasta))
    if not Path(read_sequences_fasta).exists():
        print("Read sequences file not found at {}".format(read_sequences_fasta))

    cores = multiprocessing.cpu_count()
    cores_to_use = max(1, cores-1)
    bowtie_indexes_path = Path(inter_path("genomes/{}".format(genome_fasta_path.stem)))
    os.makedirs(bowtie_indexes_path.parent.resolve(), exist_ok=True)

    # build a Bowtie2 index for mapping from reference genome
    build_command = 'bowtie2-build {} {} -q'.format(genome_fasta_path.resolve(), bowtie_indexes_path)
    subprocess.run(build_command, shell=True)

    # align reads, output SAM alignments. See bowtie2 documentation for ful details on arguments used
    # -a mode searches for all alignments
    # -p runs on multiple parallel processors to boost performance
    # -t prints time information
    # -f to specify fasta input (ignore qualities)
    output_alignment_results = inter_path("{}_bwt2_full.sam".format(output_filename))
    align_command = 'bowtie2 -x {} -t -f {} -S {} -p {} -a --quiet'.format(bowtie_indexes_path, read_sequences_fasta, output_alignment_results, cores_to_use)
    subprocess.run(align_command, shell=True)

    # Use files with the awk because windows is crazy
    # about escaping characters in shell commands
    no_matches_cmd_path = Path('./awk_commands/no_matches')
    matches_cmd_path = Path('./awk_commands/yes_matches')

    matches_output_path = inter_path("{}_bwt2_matches.sam".format(output_filename))
    no_matches_output_path = inter_path("{}_bwt2_no_matches.sam".format(output_filename))
    make_matches_command = '''awk -f {} {} > {}'''.format(matches_cmd_path, output_alignment_results, matches_output_path)
    make_no_matches_command = '''awk -f {} {} > {}'''.format(no_matches_cmd_path, output_alignment_results, no_matches_output_path)
    subprocess.run(make_matches_command, shell=True)
    subprocess.run(make_no_matches_command, shell=True)
    return

def sam_to_fasta(sam_path, fasta_path):
    with open(sam_path, 'r') as sam_file_handler:
        with open(fasta_path, 'w') as output_fasta:
            reader = csv.reader(sam_file_handler, delimiter="\t")
            writer = csv.writer(output_fasta)
            for row in reader:
                output_fasta.write('>{}\n{}\n'.format(row[0], row[9]))
    return

def run_alignment(fingerprinted_path, meta_info):
    print("Running alignment mapping...")
    start = time.perf_counter()

    genome_file_path = Path(meta_info['Genome fasta file'])
    plasmid_file_path = Path(meta_info['Plasmid fasta file'])

    # Align to reference genome based on path provided
    genome_reads = inter_path("genome_bwt2_matches.sam").format(meta_info['Sample'])
    genome_no_reads = inter_path("genome_bwt2_no_matches.sam").format(meta_info['Sample'])
    hist_results = ''
    if not Path(genome_reads).exists() or True:
        if genome_file_path.exists():
            find_alignments(fingerprinted_path, meta_info['Genome fasta file'], 'genome')
            hist_results = correct_output_reads(genome_reads, genome_no_reads, meta_info, 'genome')
        else:
            print("No genome fasta provided in the info csv, skipping genome alignment")

    elapsed_time = round(time.perf_counter() - start, 2)
    print("Genome histogram data exists ({} seconds)".format(elapsed_time))

    # Turn genome_no_reads sam file to fasta format to run against the plasmid
    genome_no_reads_fasta = inter_path("{}.fasta".format(Path(genome_no_reads).stem))
    sam_to_fasta(genome_no_reads, genome_no_reads_fasta)

    # If a path for reference plasmid was provided, align to plasmid
    if len(meta_info['Plasmid fasta file']) > 1:
        plasmid_reads = inter_path("plasmid_bwt2_matches.sam").format(meta_info['Sample'])
        plasmid_no_reads = inter_path("plasmid_bwt2_no_matches.sam").format(meta_info['Sample'])
        if not Path(plasmid_reads).exists() or True:
            if plasmid_file_path.exists():
                find_alignments(genome_no_reads_fasta, meta_info['Plasmid fasta file'], 'plasmid')
                correct_output_reads(plasmid_reads, plasmid_no_reads, meta_info, 'plasmid')
            else:
                print("No plasmid fasta provided in the csv, skipping plasmid alignment")
    
        elapsed_time = round(time.perf_counter() - start, 2)
        print("Plasmid histogram data exists ({} seconds)".format(elapsed_time))

    return hist_results

def sam_to_chunks(csv_file):  # breaks down large sam file into smaller chunks to reduce memory usage
    chunks = pd.read_csv(csv_file,sep="\t",usecols=[0,1,2,3,4,9,11,12,13,15,16,17],header=None, chunksize=1000000)
    return chunks

def correct_reads(matches_sam, output_name):
    '''
    # dictionary parsing reads
    reads = {}
    with sam_file_handler as open(matches_sam, 'r'):
        reader = csv.reader(sam_file_handler, delimited="\t")
        for row in reader:
            read_number = row[0]
            coordinate = row[3]
            flag_sum = row[2]
            if flag_sum == 0:
                direction = 'LR'
                corrected_coor = coordinate + map_length
            else:
                direction = 'RL'
                corrected_coor = coordinate + TSD
            min_obj = {'coord': corrected_coor, 'direction': direction}
            
            if reads[read_number]:
                reads.read_number.push(min_obj)
            else:
                reads = [min_obj]

    unique_reads = [read for read in reads if read.length == 1]
    non_unique_reads = [read fo]
    '''
    col_names = "read_number, flag_sum, ref_genome, ref_genome_coordinate, mapq, read_sequence, AS, XN, XM, XG, NM, MD".split(", ")
    
    SAM_full = sam_to_chunks(matches_sam)
    unique_read_numbers = []
    for chunk in SAM_full:
        chunk.columns = col_names
        chunk_unique_read_numbers = chunk.read_number.value_counts()[chunk.read_number.value_counts() == 1]
        unique_read_numbers += list(chunk_unique_read_numbers.index)

        if len(set(unique_read_numbers)) != len(unique_read_numbers):
            print("Annoying, there was a duplicate between chunks, finding and eliminating it now...")
            unique_read_numbers = [num for num in unique_read_numbers if unique_read_numbers.count(num) == 1]

    unique_reads = pd.DataFrame(columns=col_names)
    non_unique_reads = pd.DataFrame(columns=col_names)
    for chunk in sam_to_chunks(matches_sam):
        chunk.columns = col_names
        chunk_unique_reads = chunk[chunk.read_number.isin(unique_read_numbers)]
        unique_reads = unique_reads.append(chunk_unique_reads)
        chunk_non_unique_reads = chunk[np.logical_not(chunk.read_number.isin(unique_read_numbers))]
        non_unique_reads = non_unique_reads.append(chunk_non_unique_reads)

    #Get the corrected integration coordinate

    histogram = unique_reads[['read_number','flag_sum','ref_genome_coordinate','read_sequence']]

    corrected_coor = []
    orientation = []
    for i,j in zip(histogram.ref_genome_coordinate, histogram.flag_sum):
        if j == 0:
            # j=0 means the FP was matched on the forward strand
            corrected_coor.append(i + map_length)
            orientation.append('LR')
        else:
            # this means it was matched on the reverse complement
            corrected_coor.append(i + TSD)
            orientation.append('RL')

    histogram = histogram.assign(corrected_coor=corrected_coor, orientation=orientation)
    
    RL_counts = histogram[histogram.orientation == 'RL'].corrected_coor.value_counts().sort_index(0)
    LR_counts = histogram[histogram.orientation == 'LR'].corrected_coor.value_counts().sort_index(0)
    # print("RL", RL_counts)
    # print("LR", LR_counts)

    counts = histogram.corrected_coor.value_counts().sort_index(0)
    counts.index.name = 'position'
    counts.name = 'reads'
    # Decrement to make the counts all the exact same as what Biopython sequence.seq.find() would give
    # Presumably because bowtie is 1-indexing and biopython is 0-indexing
    counts.index -= 1
    
    hist_path = output_path("{}_read_locations.csv".format(output_name))
    counts.to_csv(hist_path)
    return [histogram, unique_reads, non_unique_reads]

def correct_output_reads(matches_sam, no_matches_sam, meta_info, output_name):
    donor_fp = meta_info['Donor sequence']
    spike_fp = meta_info['Spike in sequence']

    histogram, unique_reads, non_unique_reads = correct_reads(matches_sam, output_name)
    
    fasta_sequence = []
    for i,j,k in zip(histogram.read_number,histogram.read_sequence,histogram.flag_sum):
        if k !=0:
            fasta_sequence.append(">%s"%(i) + "\n" + str(Seq(j).reverse_complement()))
        else:
            fasta_sequence.append(">%s"%(i) + "\n" + j)

    fasta_file = "\n".join(fasta_sequence)
    fasta_path = output_path("{}_unique_reads.fasta".format(output_name))
    with open(fasta_path, 'w', newline='') as file:
        file.write(fasta_file)

    ## Check the fingerprints without genome matches against donor, spike and CRISPR Sequence
    ##
    ## We check the sequences with no genome matches against these donor and spike first. 
    ## If it doesn't match either of them, we additionally check a
    ## sample-specific CRISPR Sequence Array,
    ## These read counts are added to the output logs for the run. 
    no_match_sequences = pd.read_csv(no_matches_sam, sep="\t", usecols=[0,1,2,3,4,9,11,12,13,15,16,17], header=None)
    col_names = "read_number, flag_sum, ref_genome, ref_genome_coordinate, mapq, read_sequence, AS, XN, XM, XG, NM, MD".split(", ")
    no_match_sequences.columns = col_names

    crispr_array_seq = Seq(meta_info['CRISPR Array Sequence']).upper()
    crispr_array_seq_rc = crispr_array_seq.reverse_complement()
    donor_matches = 0
    spike_matches = 0
    cripsr_seq_matches = 0
    for read_seq in no_match_sequences['read_sequence']:
        if read_seq == donor_fp:
            donor_matches += 1
        elif read_seq == spike_fp:
            spike_matches += 1
        elif crispr_array_seq and (read_seq in crispr_array_seq or read_seq in crispr_array_seq_rc):
            cripsr_seq_matches += 1

    unique_reads_seq_count = len(unique_reads.read_number.unique())
    non_unique_reads_seq_count = len(non_unique_reads.read_number.unique())
    output = {
        'Unique Genome-Mapping Reads': unique_reads_seq_count,
        'Total Genome-Mapping Reads': non_unique_reads_seq_count + unique_reads_seq_count,
        'Undigested Donor Reads': donor_matches,
        'Spike-in Reads': spike_matches,
        'CRISPR Array Self-Targeting Reads': cripsr_seq_matches
    }
    print("Alignment Summary------------------")
    for stat in output:
        print(stat,':', output[stat])
    print("-----------------------------------")
    return output