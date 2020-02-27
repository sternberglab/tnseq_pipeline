from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import time
import os

from parameters import flank_sequence, fingerprint_length, delete_intermediates

def hamming_dist(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def fpgen(input_reads):  # generator function that returns
    flank_length = len(flank_sequence)
    for record in SeqIO.parse(input_reads, "fastq"):
        i = 0
        while 0 <= i <= (len(record.seq) - flank_length):
            query = record.seq[i:i+flank_length]
            if hamming_dist(query, flank_sequence) < 2:
                if i >= fingerprint_length:
                    fingerprint = record.seq[i - fingerprint_length:i]
                    newrec = SeqRecord(fingerprint, id=record.id, name=record.name)
                    yield newrec
                i += 2000
            else:
                i += 1

def fingerprinting(input_reads, fingerprinted_path):  # main command function for filtering
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
