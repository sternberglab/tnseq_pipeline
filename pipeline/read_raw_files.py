from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time
import os

from .utils import inter_path
from parameters import Qscore_threshold, delete_intermediates

def filtergen(file, threshold):  # generator function that returns edited reads that pass filter, to write new fastq file
    for record in SeqIO.parse(file, "fastq"):
        # Convert base qualities to Boolean based on Qscore threshold value. Only use reads with >=50% non-N:
        recordqual = [x > threshold for x in record.letter_annotations['phred_quality']]  # list of True, False etc
        if float(sum(recordqual)) / float(len(recordqual)) >= .5:  # note that True = 1, False = 0 for summing
            # create new SeqRecord with edited read sequence
            newrec = SeqRecord(record.seq, id=record.id, name=record.name,
                               description=record.description, letter_annotations=record.letter_annotations)
            yield newrec

def process_files(code, input_filenames, filtered_path):
    print("Processing and quality filtering {} files...".format(len(input_filenames)))
    start = time.perf_counter()
    total_records = 0
    filtered_records = 0
    concat_path = inter_path('{}_Q{}_CONCAT.fastq'.format(code, Qscore_threshold))
    with open(inter_path(concat_path), 'w') as concat_outfile:
        for file in input_filenames:
            with open(file, 'r') as fastq_file:
                for line in fastq_file:
                    concat_outfile.write(line)
                    if line.startswith('@'):
                        total_records += 1
        
    with open(filtered_path, 'w+') as filtered_outfile:
        filtered_records = SeqIO.write(filtergen(concat_path, Qscore_threshold), filtered_outfile, "fastq")
    
    if delete_intermediates:
        os.remove(concat_path)

    elapsed = round(time.perf_counter() - start, 2)
    percent_passing = round((total_records - filtered_records) / total_records, 2)
    print("Read {} total records, {} filtered records ({}%) in {} seconds".format(total_records, filtered_records, percent_passing, elapsed))
    return {'raw_record_count': total_records, 'filtered_record_count': filtered_records}