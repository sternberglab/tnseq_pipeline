from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time
import os
import zipfile
from pathlib import Path

from .utils import inter_path
from parameters import Qscore_threshold, delete_intermediates, working_dir

total_records = 0
def filtergen(files, threshold):  # generator function that returns edited reads that pass filter, to write new fastq file
    global total_records
    total_records = 0
    for file in files:
        for record in SeqIO.parse(file, "fastq"):
            total_records += 1
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
    passing_filter = 0
    with open(filtered_path, 'w+') as filtered_outfile:
        passing_filter = SeqIO.write(filtergen(input_filenames, Qscore_threshold), filtered_outfile, "fastq")

    elapsed = round(time.perf_counter() - start, 2)
    percent_passing = 100 - round(100 * (total_records - passing_filter) / total_records, 2)
    print("Read {} total records, {} filtered records ({}%) in {} seconds".format(total_records, passing_filter, percent_passing, elapsed))
    return {'Total Raw Reads': total_records, 'Filtered Reads': passing_filter}

def unzip_files():
    raw_files_dir = os.path.join(Path(working_dir), 'raw')
    zip_files = Path(raw_files_dir).glob('*.zip')
    for zipped in zip_files:
        with zipfile.ZipFile(zipped.resolve(), 'r') as zip_ref:
            zip_ref.extractall(raw_files_dir)
        if delete_intermediates:
            os.remove(zipped)