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
from parameters import transposon_end_flanking_sequence_length as map_length, target_site_duplication_length as TSD

def find_alignments(read_sequences_fasta, target_fasta, output_filename):
    print("Find alignments for {}".format(output_filename))

    target_fasta_path = Path(target_fasta)
    if not target_fasta_path.exists():
        print("Source sequence file not found at {}".format(target_fasta))
    if not Path(read_sequences_fasta).exists():
        print("Read sequences file not found at {}".format(read_sequences_fasta))

    cores = multiprocessing.cpu_count()
    cores_to_use = max(1, cores-1)
    bowtie_indexes_path = Path(inter_path("genomes/{}".format(target_fasta_path.stem))).resolve()
    os.makedirs(bowtie_indexes_path.parent.resolve(), exist_ok=True)
    build_command = 'bowtie2-build "{}" "{}"'.format(target_fasta_path.resolve(), bowtie_indexes_path)
    subprocess.run(build_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    output_alignment_results = inter_path("{}_bwt2_full.sam".format(output_filename))
    align_command = 'bowtie2 -x {} -f {} -S {} -p {} -a --quiet'.format(bowtie_indexes_path, read_sequences_fasta, output_alignment_results, cores_to_use)
    subprocess.run(align_command, shell=True)

    matches_output_path = inter_path("{}_bwt2_matches.sam".format(output_filename))
    no_matches_output_path = inter_path("{}_bwt2_no_matches.sam".format(output_filename))
    matches_output_file = open(matches_output_path, 'w')
    no_matches_output_file = open(no_matches_output_path, 'w')
    with open(output_alignment_results, 'r') as full_sam:
        for line in full_sam:
            if line[0] == '@':
                continue
            if 'NM:i:0' in line:
                matches_output_file.write(line)
            else:
                no_matches_output_file.write(line)

    matches_output_file.close()
    no_matches_output_file.close()
    return

def sam_to_fasta(sam_path, fasta_path):
    with open(sam_path, 'r') as sam_file_handler:
        with open(fasta_path, 'w') as output_fasta:
            reader = csv.reader(sam_file_handler, delimiter="\t")
            writer = csv.writer(output_fasta)
            for row in reader:
                output_fasta.write('>{}\n{}\n'.format(row[0], row[9]))
    return

def run_alignment(flanks_path, meta_info):
    print("Running alignment mapping...")
    start = time.perf_counter()

    target_file = meta_info['Target fasta file']
    target_file_path = Path(target_file)
    second_target_file = meta_info['Second target fasta file']
    second_target_file_path = Path(second_target_file)

    target_reads = inter_path(f"{meta_info['Sample']}_target_bwt2_matches.sam")
    target_no_reads = inter_path(f"{meta_info['Sample']}_target_bwt2_no_matches.sam")
    hist_results = ''
    if len(target_file) > 1:
        find_alignments(flanks_path, target_file_path, f"{meta_info['Sample']}_target")
        hist_results = correct_output_reads(target_reads, target_no_reads, meta_info, f"{meta_info['Sample']}_target")
    else:
        print("No target fasta provided in the info csv, skipping genome alignment")

    elapsed_time = round(time.perf_counter() - start, 2)
    print("Target data processed in ({} seconds)".format(elapsed_time))
    # Turn target_no_reads sam file to fasta format to run against the second target
    target_no_reads_fasta = inter_path("{}.fasta".format(Path(target_no_reads).stem))
    sam_to_fasta(target_no_reads, target_no_reads_fasta)

    if len(second_target_file) > 1:
        second_target_reads = inter_path(f"{meta_info['Sample']}_second_target_bwt2_matches.sam")
        second_target_no_reads = inter_path(f"{meta_info['Sample']}_second_target_bwt2_no_matches.sam")
        if second_target_file_path.exists():
            find_alignments(target_no_reads_fasta, second_target_file_path, f"{meta_info['Sample']}_second_target")
            second_results = correct_output_reads(second_target_reads, second_target_no_reads, meta_info, f"{meta_info['Sample']}_second_target")
            hist_results.update(second_results)
        else:
            print("No second target fasta provided in the csv, skipping second target alignment")
    
        elapsed_time = round(time.perf_counter() - start, 2)
        print("Second target data processed in ({} seconds)".format(elapsed_time))

    return hist_results

def sam_to_chunks(csv_file):
    chunks = pd.read_csv(csv_file,sep="\t",usecols=[0,1,2,3,4,9,11,12,13,15,16,17],header=None, chunksize=1000000)
    return chunks


def correct_read(genome_coord, read_is_fw_strand, spacer_is_fw_strand, corrected_coor, orientation, spacer_coord, default_orientation):
    alt_orientation = 'LR' if default_orientation == 'RL' else 'RL'
    if read_is_fw_strand:
        # j=0 means the FP was matched on the forward strand
        # j=256 means still forward strand, was a secondary alignment
        # READ HERE: CODE THAT ASSIGNS RL/LR AND T'RL/T'LR AND GENOMIC COORDINATES
        if not spacer_coord:
            coord = genome_coord + map_length
            read_orient = default_orientation
        elif not spacer_is_fw_strand:
            if spacer_coord > (
                    genome_coord + map_length - TSD):  # if T-LR and spacer on reverse strand, read on FW strand, need spacer coord to be greater than read coord
                coord = genome_coord + map_length - TSD
            else:
                coord = genome_coord + map_length
            read_orient = alt_orientation
        else:
            if spacer_coord < (genome_coord + map_length):
                coord = genome_coord + map_length
            else:
                coord = genome_coord + map_length - TSD
            read_orient = default_orientation
        corrected_coor.append(coord)
        orientation.append(read_orient)
    
    # this means it was matched on the reverse complement
    else:
        if not spacer_coord:
            coord = genome_coord
            read_orient = default_orientation
        elif not spacer_is_fw_strand:
            if spacer_coord and spacer_coord > genome_coord:
                coord = genome_coord
            else:
                coord = genome_coord + TSD
            read_orient = default_orientation
        else:
            if spacer_coord and spacer_coord < genome_coord:
                coord = genome_coord + TSD
            else:
                coord = genome_coord
            read_orient = alt_orientation
        corrected_coor.append(coord)
        orientation.append(read_orient)

def correct_reads(matches_sam, output_name, meta_info):
    if 'second' not in output_name:
        genome = SeqIO.read(Path(meta_info['Target fasta file']), "fasta")
    else:
        genome = SeqIO.read(Path(meta_info['Second target fasta file']), "fasta")
    refseq = genome.seq.upper()
    rev_refseq = genome.seq.reverse_complement().upper()
    spacer_seq = meta_info['Spacer'].upper()
    spacer_is_fw_strand = refseq.find(spacer_seq) >= 0 or rev_refseq.find(spacer_seq) == -1
    if spacer_is_fw_strand:
        if refseq.find(spacer_seq) >= 0:
            spacer_end_coord = int(refseq.find(spacer_seq)) + len(spacer_seq)
            meta_info['SpacerStartRefSeq'] = f"fw_{refseq.find(spacer_seq)}"
        else: 
            spacer_end_coord = None
            meta_info['SpacerStartRefSeq'] = None
    else:
        if rev_refseq.find(spacer_seq) >= 0:
            spacer_end_coord = len(genome.seq) - rev_refseq.find(spacer_seq) - len(spacer_seq)
            meta_info['SpacerStartRefSeq'] = f"rv_{len(genome.seq) - rev_refseq.find(spacer_seq)}"
        else: 
            spacer_end_coord = None
            meta_info['SpacerStartRefSeq'] = None
    
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
        unique_reads = pd.concat([unique_reads, chunk_unique_reads])
        chunk_non_unique_reads = chunk[np.logical_not(chunk.read_number.isin(unique_read_numbers))]
        non_unique_reads = pd.concat([non_unique_reads, chunk_non_unique_reads])

    #Get the corrected integration coordinate
    histogram = unique_reads[['read_number','flag_sum','ref_genome_coordinate','read_sequence']]

    is_left_end = meta_info.get('transposon_end_side', 'Right').lower() == 'left'
    default_orientation = 'LR' if is_left_end else 'RL'
    corrected_coor = []
    orientation = []
    n=0
    k=0
    for i,j,m in zip(histogram.ref_genome_coordinate, histogram.flag_sum, histogram.read_sequence):
        read_is_fw_strand = j%256 == 0
        if meta_info['read_type'] == 'fragment':
            # fragmentation reads use the RC of the tn_end sequence, 
            # so their "actual" tn_end sequence read would be the opposite
            # strand of the detected one
            read_is_fw_strand = not read_is_fw_strand
        correct_read(i, read_is_fw_strand, spacer_is_fw_strand, corrected_coor, orientation, spacer_end_coord, default_orientation):
    histogram = histogram.assign(corrected_coor=corrected_coor, orientation=orientation)

    RL_counts = histogram[histogram.orientation == 'RL'].corrected_coor.value_counts().sort_index()
    LR_counts = histogram[histogram.orientation == 'LR'].corrected_coor.value_counts().sort_index()
    totals = histogram.corrected_coor.value_counts().sort_index()

    RL_counts.name = "RL"
    LR_counts.name = "LR"
    totals.name = 'reads'

    # Decrement to make the counts all the exact same as what Biopython sequence.seq.find() would give
    # Presumably because bowtie is 1-indexing and biopython is 0-indexing
    combine = pd.concat([totals, RL_counts, LR_counts], axis=1)
    combine.index.name = "position"
    combine.index -= 1
    combine["RL"] = combine["RL"].fillna(0)
    combine["LR"] = combine["LR"].fillna(0)
    
    hist_path = output_path(os.path.join('samples', f"{meta_info['output_date']}_{output_name}_read_locations.csv"))
    combine.to_csv(hist_path)
    return [histogram, unique_reads, non_unique_reads]

def correct_output_reads(matches_sam, no_matches_sam, meta_info, output_name):
    donor_fp = meta_info['Donor sequence']
    spike_fp = meta_info['Spike in sequence']

    try:
        histogram, unique_reads, non_unique_reads = correct_reads(matches_sam, output_name, meta_info)

        fasta_sequence = []
        for i,j,k in zip(histogram.read_number,histogram.read_sequence,histogram.flag_sum):
            if k !=0:
                fasta_sequence.append(">%s"%(i) + "\n" + str(Seq(j).reverse_complement()))
            else:
                fasta_sequence.append(">%s"%(i) + "\n" + j)

        fasta_file = "\n".join(fasta_sequence)
        fasta_path = inter_path("{}_unique_reads.fasta".format(output_name))
        with open(fasta_path, 'w', newline='') as file:
            file.write(fasta_file)
        unique_reads_seq_count = len(unique_reads.read_number.unique())
        non_unique_reads_seq_count = len(non_unique_reads.read_number.unique())
    except:
        unique_reads_seq_count = 0
        non_unique_reads_seq_count = 0
    
    ## Check the flanking sequences without genome matches against donor, spike and CRISPR Sequence
    ##
    ## We check the sequences with no genome matches against these donor and spike first. 
    ## If it doesn't match either of them, we additionally check a
    ## sample-specific CRISPR Sequence Array,
    ## These read counts are added to the output logs for the run. 
    no_match_sequences = pd.read_csv(no_matches_sam, sep="\t", usecols=[0,1,2,3,4,9], header=None)
    col_names = "read_number, flag_sum, ref_genome, ref_genome_coordinate, mapq, read_sequence".split(", ")
    no_match_sequences.columns = col_names

    crispr_array_seq = meta_info.get('CRISPR Array Sequence', None)
    crispr_array_seq_rc = None
    if crispr_array_seq:
        crispr_array_seq = Seq(meta_info['CRISPR Array Sequence']).upper() if 'CRISPR Array Sequence' in meta_info else None
        crispr_array_seq_rc = crispr_array_seq.reverse_complement() if crispr_array_seq else None
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

    output = {}
    if 'second' not in output_name:
        output = {
            'Unique Target Mapping Reads': unique_reads_seq_count,
            'Total Target Mapping Reads': non_unique_reads_seq_count + unique_reads_seq_count,
            'Contaminating donor reads': donor_matches,
            'Spike-in Reads': spike_matches,
            'CRISPR Array Self-Targeting Reads': cripsr_seq_matches
        }
    else:
        output = {
            'Unique Second Target Mapping Reads': unique_reads_seq_count,
            'Total Second Target Mapping Reads': non_unique_reads_seq_count + unique_reads_seq_count
        }
    return output
