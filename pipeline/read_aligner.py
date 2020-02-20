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

from .utils import inter_path
from parameters import delete_intermediates, genome_path, fingerprint_length as map_length, transposon_site_duplication_length as TSD

def run_alignment(fingerprinted_path, run_prefix):
    print("Running alignment mapping...")
    output_no_mismatch_sam = inter_path("{}_bwt2_no_mismatches_full.sam".format(run_prefix))
    
    if not Path(output_no_mismatch_sam).exists():
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
        
        filter_command = '''cat {} | awk '$0 ~"NM:i:0"' > {}'''.format(output_full_sam, output_no_mismatch_sam)
        subprocess.run(filter_command, shell=True)
        
        if delete_intermediates:
            shutil.rmtree(bowtie_indexes_path.parent.resolve())
            os.remove(output_full_sam)

    print("Generating the histogram data...")
    hist_results = make_histogram(output_no_mismatch_sam, run_prefix)
    return hist_results

def make_histogram(sam_path, run_prefix):
    SAM_full = pd.read_csv(sam_path,sep="\t",usecols=[0,1,2,3,4,9,11,12,13,15,16,17],header=None)
    col_names = "read_number, flag_sum, ref_genome, ref_genome_coordinate, mapq, read_sequence, AS, XN, XM, XG, NM, MD".split(", ")
    SAM_full.columns = col_names

    unique_read_numbers = SAM_full.read_number.value_counts()[SAM_full.read_number.value_counts() == 1]
    unique_reads = SAM_full[SAM_full.read_number.isin(list(unique_read_numbers.index))]
    non_unique_reads = SAM_full[np.logical_not(SAM_full.read_number.isin(unique_read_numbers.index))]

    #Get the corrected integration coordinate
    histogram = unique_reads[['read_number','flag_sum','ref_genome_coordinate','read_sequence']]

    corrected_coor = []
    for i,j in zip(histogram.ref_genome_coordinate,histogram.flag_sum):
        if j == 0:
            corrected_coor.append(i + map_length)
        else:
            corrected_coor.append(i + TSD)

    histogram['corrected_coor'] = corrected_coor
    counts = histogram.corrected_coor.value_counts()

    histogram_count = []
    genome_length = len(SeqIO.read(Path(genome_path).resolve(), 'fasta').seq)
    for i in range(genome_length):
        if i in counts:
            histogram_count.append(counts[i])
        else:
            histogram_count.append(0)

    hist = pd.DataFrame(histogram_count,columns=['count'])
    hist.index = range(1,len(hist)+1)
    hist.index.name = 'position'
    hist_path = inter_path("{}_unique_reads_aligned_histogram.csv".format(run_prefix))
    hist.to_csv(hist_path)

    fasta_sequence = []
    for i,j,k in zip(histogram.read_number,histogram.read_sequence,histogram.flag_sum):
        if k !=0:
            fasta_sequence.append(">%s"%(i) + "\n" + str(Seq(j).reverse_complement()))
        else:
            fasta_sequence.append(">%s"%(i) + "\n" + j)

    fasta_file = "\n".join(fasta_sequence)
    fasta_path = inter_path("{}_unique_reads.fasta".format(run_prefix))
    with open(fasta_path, 'w', newline='') as file:
        file.write(fasta_file)
    return {
        'unique_reads_count': len(unique_reads.read_number.unique()),
        'non_unique_reads_count': len(non_unique_reads.read_number.unique())
    }
