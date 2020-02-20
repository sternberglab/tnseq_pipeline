from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path

from parameters import flank_sequence, fingerprint_length

def hamming_dist(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def fpgen(filename):  # generator function that returns
    flank_length = len(flank_sequence)
    for record in SeqIO.parse(filename, "fastq"):
        i = 0
        while 0 <= i <= (len(record.seq) - flank_length):
            query = record.seq[i:i+flank_length]
            if hamming_dist(query, flank_sequence) < 2:
                if i >= fp_length:
                    fingerprint = record.seq[i - fingerprint_length:i]
                    newrec = SeqRecord(fingerprint, id=record.id, name=record.name)
                    yield newrec
                i += 2000
            else:
                i += 1

def fingerprinting(filename, outfile):  # main command function for filtering
    total = 0
    # generate input records to write new fingerprinted fastq
    with open(outfile, 'w+') as fastafile:
        fp_total = SeqIO.write(fpgen(filename), fastafile, "fasta")

    return (total, fp_total)
