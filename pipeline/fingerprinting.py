from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import time
import os

from parameters import flank_sequence, fingerprint_length, delete_intermediates

# Function that determines Hamming distance between two strings
# Useful for searching for reads with mismatches
def hamming_dist(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

# Generator function that yields a fingerprint sequence for each read
# Uses the transposon end sequence (flank_sequence) to locate transposon within read
# Then checks whether FP sequence upstream is of sufficient length, then yields the FP sequence
# FP sequence records do not contain PHRED quality scores; output to fasta only.
def fpgen(input_reads):
    flank_length = len(flank_sequence)
    for record in SeqIO.parse(input_reads, "fastq"):
        found = record.seq.find(flank_sequence)
        # First search for perfect match for flank_sequence within read
        # If not found, uses Hamming dist to search while allowing mismatches
        if found > -1:
            if found >= fingerprint_length:
                fingerprint = record.seq[found - fingerprint_length:found]
                newrec = SeqRecord(fingerprint, id=record.id, name=record.name)
                yield newrec
        
        else:  # iterate through each window of the read and determines hamm dist
            i = 0
            while 0 <= i <= (len(record.seq) - flank_length):
                query = record.seq[i:i+flank_length]
                if hamming_dist(query, flank_sequence) < 2:  # matching with one mismatch
                    if i >= fingerprint_length:
                        fingerprint = record.seq[i - fingerprint_length:i]
                        newrec = SeqRecord(fingerprint, id=record.id, name=record.name)
                        yield newrec
                    i += 2000
                else:
                    i += 1

# Function that calls fpgen to output fingerprints into a fasta file
def fingerprinting(input_reads, fingerprinted_path):
    print("Fingerprinting records...")
    start = time.perf_counter()
    # generate input records to write new fingerprinted fastq
    with open(fingerprinted_path, 'w+') as fastafile:
        fp_reads = SeqIO.write(fpgen(input_reads), fastafile, "fasta")
    
    elapsed = round(time.perf_counter() - start, 2)
    print("{} reads had a valid fingerprint (ran in {} seconds)".format(fp_reads, elapsed))
    if delete_intermediates:
        os.remove(input_reads)
    return {'Valid Fingerprint Reads': fp_reads}
