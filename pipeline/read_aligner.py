import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import pandas as pd
import sys
import multiprocessing
import shutil
import numpy as np
import subprocess
import platform
import time

from .utils import inter_path, output_path
from parameters import delete_intermediates, genome_path, fingerprint_length as map_length, transposon_site_duplication_length as TSD

# Unique parameters for this file
donor_fp = "TGCTGAAACCTCAGGCA"
spike_fp = "ATGATTACGCCAAGCTT"


def run_alignment(fingerprinted_path, run_prefix, meta_info):
    print("Running alignment mapping...")
    start = time.perf_counter()
    output_no_mismatch_sam = inter_path("{}_bwt2_no_mismatches_full.sam".format(run_prefix))
    output_mismatch_sam = inter_path("{}_bwt2_mismatches_full.sam".format(run_prefix))

    if not Path(output_no_mismatch_sam).exists() or True:
        print("Finding genome alignments...")
        genome_file = Path(genome_path)
        # First have bowtie build and index the genome if it hasn't yet
        if not genome_file.exists():
            raise ValueError("Genome file can't be found at {}".format(file))

        cores = multiprocessing.cpu_count()
        cores_to_use = max(1, cores-1)

        bowtie_indexes_path = Path(inter_path("genomes/{}".format(genome_file.stem)))
        os.makedirs(bowtie_indexes_path.parent.resolve(), exist_ok=True)
        
        build_command = 'bowtie2-build {} {} -q'.format(genome_file.resolve(), bowtie_indexes_path)
        subprocess.run(build_command, shell=True)

        output_full_sam = inter_path("{}_bwt2_full.sam".format(run_prefix))
        align_command = 'bowtie2 -x {} -t -f {} -S {} -p {} -a --quiet'.format(bowtie_indexes_path, fingerprinted_path, output_full_sam, cores_to_use)
        subprocess.run(align_command, shell=True)
        
        # use "cat" for non-windows, "type" for windows
        if platform.system() == 'Windows':
            make_matches_command = '''awk "$0 ~\""NM:i:0\""" {} > {}'''.format(output_full_sam, output_no_mismatch_sam)
            make_no_matches_command = '''awk "$0 !~\""NM:i:0\"" && $0 !~ /^@/" {} > {}'''.format(output_full_sam, output_mismatch_sam)
        else:
            make_matches_command = '''awk '$0 ~"NM:i:0"' {} > {}'''.format(output_full_sam, output_no_mismatch_sam)
            make_no_matches_command = '''awk '$0 !~"NM:i:0" && $0 !~ /^@/' {} > {}'''.format(output_full_sam, output_mismatch_sam)

        subprocess.run(make_matches_command, shell=True)
        subprocess.run(make_no_matches_command, shell=True)       

    print("Generating the histogram data...")
    hist_results = correct_output_reads(output_no_mismatch_sam, output_mismatch_sam, meta_info, run_prefix)
    elapsed_time = round(time.perf_counter() - start, 2)
    print("Finished doing genome mapping and generating histogram data in {} seconds".format(elapsed_time))
    if delete_intermediates:
        os.remove(output_mismatch_sam)
        os.remove(output_no_mismatch_sam)
        shutil.rmtree(bowtie_indexes_path.parent.resolve())
        os.remove(output_full_sam)
    return hist_results

def correct_output_reads(matches_sam, no_matches_sam, meta_info, run_prefix):
    SAM_full = pd.read_csv(matches_sam,sep="\t",usecols=[0,1,2,3,4,9,11,12,13,15,16,17],header=None)
    col_names = "read_number, flag_sum, ref_genome, ref_genome_coordinate, mapq, read_sequence, AS, XN, XM, XG, NM, MD".split(", ")
    SAM_full.columns = col_names

    unique_read_numbers = SAM_full.read_number.value_counts()[SAM_full.read_number.value_counts() == 1]
    unique_reads = SAM_full[SAM_full.read_number.isin(list(unique_read_numbers.index))]
    non_unique_reads = SAM_full[np.logical_not(SAM_full.read_number.isin(unique_read_numbers.index))]

    #Get the corrected integration coordinate
    histogram = unique_reads[['read_number','flag_sum','ref_genome_coordinate','read_sequence']]

    corrected_coor = []
    for i,j in zip(histogram.ref_genome_coordinate, histogram.flag_sum):
        if j == 0:
            corrected_coor.append(i + map_length)
        else:
            corrected_coor.append(i + TSD)

    histogram['corrected_coor'] = corrected_coor

    counts = histogram.corrected_coor.value_counts().sort_index(0)
    counts.index.name = 'position'
    # Increment to make the counts all the exact same as original analysis
    # Presumably because the old script was 1-indexing and bowtie is 0-indexing
    counts.index += 1
    
    hist_path = output_path("{}_genome_read_locations.csv".format(run_prefix))
    counts.to_csv(hist_path)
    
    genome_length = len(SeqIO.read(Path(genome_path).resolve(), 'fasta').seq)

    fasta_sequence = []
    for i,j,k in zip(histogram.read_number,histogram.read_sequence,histogram.flag_sum):
        if k !=0:
            fasta_sequence.append(">%s"%(i) + "\n" + str(Seq(j).reverse_complement()))
        else:
            fasta_sequence.append(">%s"%(i) + "\n" + j)

    fasta_file = "\n".join(fasta_sequence)
    fasta_path = output_path("{}_unique_reads.fasta".format(run_prefix))
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

    crispr_array_seq = meta_info['CRISPR Array Sequence']
    donor_matches = 0
    spike_matches = 0
    cripsr_seq_matches = 0
    for read_seq in no_match_sequences['read_sequence']:
        if read_seq == donor_fp:
            donor_matches += 1
        elif read_seq == spike_fp:
            spike_matches += 1
        elif crispr_array_seq and read_seq in crispr_array_seq:
            cripsr_seq_matches += 1

    return {
        'Unique Genome-Mapping Reads': len(unique_reads.read_number.unique()),
        'Non-Unique Genome-Mapping Reads': len(non_unique_reads.read_number.unique()),
        'Donor Reads': donor_matches,
        'Spike Reads': spike_matches,
        'CRISPR Array Reads': cripsr_seq_matches
    }