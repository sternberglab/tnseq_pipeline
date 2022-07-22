from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import time
import os

from parameters import transposon_end_flanking_sequence_length

allow_tn_end_mismatch = False
def hamming_dist(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def fpgen(input_reads, meta_info):  # generator function that returns
    read_type = meta_info["read_type"]
    if read_type not in ["fragment", "restriction"]:
        raise "Read type must be 'fragment' or 'restriction'"

    transposon_end_sequence = meta_info["transposon_end_sequence"].upper()

    tn_end_length = len(transposon_end_sequence)
    for record in SeqIO.parse(input_reads, "fastq"):
        found = record.seq.find(transposon_end_sequence)
        if found > -1:
            if read_type == 'restriction':
                if found < 21 and found >= transposon_end_flanking_sequence_length:
                # restriction reads for forward, found before the R tn_end sequence
                    transposon_end_flanking_sequence = record.seq[found - transposon_end_flanking_sequence_length:found]
                    newrec = SeqRecord(transposon_end_flanking_sequence, id=record.id, name=record.name)
                    yield newrec
            else:
                # fragment reads are always coming from the inside (regardless of R or L tn_end)
                # so the transposon_end_flanking_sequence comes AFTER the tn_end
                if len(record) - (found + len(transposon_end_sequence)) >= transposon_end_flanking_sequence_length:
                    start = found + len(transposon_end_sequence)
                    end = found + len(transposon_end_sequence) + transposon_end_flanking_sequence_length
                    transposon_end_flanking_sequence = record.seq[start:end]
                    newrec = SeqRecord(transposon_end_flanking_sequence, id=record.id, name=record.name)
                    yield newrec
        
        # only allow mismatch of 1 for restriction enzyme flanking sequences
        elif read_type == 'restriction' or allow_tn_end_mismatch:
            i = 0
            while 0 <= i <= (len(record.seq) - tn_end_length):
                query = record.seq[i:i+tn_end_length]
                if hamming_dist(query, transposon_end_sequence) < 2:
                    if i >= transposon_end_flanking_sequence_length and i < 21:
                        transposon_end_flanking_sequence = record.seq[i - transposon_end_flanking_sequence_length:i]
                        newrec = SeqRecord(transposon_end_flanking_sequence, id=record.id, name=record.name)
                        yield newrec
                    i += 2000
                else:
                    i += 1

def find_flanking_sequences(input_reads, flanks_path, meta_info):  # main command function for filtering
    print("Finding transposon end flanking sequences...")
    start = time.perf_counter()
    # generate input records to write new flanking sequences fastq
    with open(flanks_path, 'w+') as fastafile:
        fp_reads = SeqIO.write(fpgen(input_reads, meta_info), fastafile, "fasta")
    
    elapsed = round(time.perf_counter() - start, 2)
    print("{} reads had a valid flanksing sequence (ran in {} seconds)".format(fp_reads, elapsed))
    return {'Transposon end-containing reads': fp_reads}
