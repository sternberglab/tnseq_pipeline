from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import time
import os

from parameters import fingerprint_length

def hamming_dist(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def fpgen(input_reads, meta_info):  # generator function that returns
    read_type = meta_info["read_type"]
    if read_type not in ["fragment", "restriction"]:
        raise "Read type must be 'fragment' or 'restriction'"

    flank_sequence = meta_info["flank_sequence"]

    flank_length = len(flank_sequence)
    for record in SeqIO.parse(input_reads, "fastq"):
        found = record.seq.find(flank_sequence)
        if found > -1:
            if read_type == 'restriction':
                # restriction reads for forward, found before the R flank sequence
                if found >= fingerprint_length:
                    fingerprint = record.seq[found - fingerprint_length:found]
                    newrec = SeqRecord(fingerprint, id=record.id, name=record.name)
                    yield newrec
            else:
                # fragment reads are always coming from the inside (regardless of R or L flank)
                # so the fingerprint comes AFTER the flank
                if len(record) - (found + len(flank_sequence)) >= fingerprint_length:
                    start = found + len(flank_sequence)
                    end = found + len(flank_sequence) + fingerprint_length
                    fingerprint = record.seq[start:end]
                    newrec = SeqRecord(fingerprint, id=record.id, name=record.name)
                    yield newrec
        
        # only allow mismatch of 1 for restriction enzyme fingerprints
        elif read_type == 'restriction':
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

def fingerprinting(input_reads, fingerprinted_path, meta_info):  # main command function for filtering
    print("Fingerprinting records...")
    start = time.perf_counter()
    # generate input records to write new fingerprinted fastq
    with open(fingerprinted_path, 'w+') as fastafile:
        fp_reads = SeqIO.write(fpgen(input_reads, meta_info), fastafile, "fasta")
    
    elapsed = round(time.perf_counter() - start, 2)
    print("{} reads had a valid fingerprint (ran in {} seconds)".format(fp_reads, elapsed))
    return {'Valid Fingerprint Reads': fp_reads}
