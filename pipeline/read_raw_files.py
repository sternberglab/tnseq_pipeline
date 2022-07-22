from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time
import gzip
import shutil
from itertools import chain
import os
import zipfile
from pathlib import Path

from .utils import inter_path
from parameters import Qscore_threshold, delete_intermediates, read_files_dir

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

def process_files(code, input_filenames, filtered_path, meta_info):
    print("Processing and quality filtering {} files...".format(len(input_filenames)))
    start = time.perf_counter()
    passing_filter = 0

    # find the read_end (R1 or R2) with the tn flank_seq and only look at it
    read_end = None
    flank_seq = meta_info['transposon_end_sequence'].upper()
    r1 = next((f for f in input_filenames if 'R1' in f.name), None)
    r1ct = 0
    if r1:
        with open(r1, 'r') as f_in:
            r1ct = f_in.read().count(flank_seq)
    r2 = next((f for f in input_filenames if 'R2' in f.name), None)
    if r2:
        with open(r2, 'r') as f_in:
            r2ct = f_in.read().count(flank_seq)
    if r1 and r2:
        if r1ct > r2ct:
            read_end = 'R1'
        else:
            read_end = 'R2'

    if read_end:
        input_filenames = [f for f in input_filenames if read_end in f.name]
    
    with open(filtered_path, 'w+') as filtered_outfile:
        passing_filter = SeqIO.write(filtergen(input_filenames, Qscore_threshold), filtered_outfile, "fastq")

    elapsed = round(time.perf_counter() - start, 2)
    percent_passing = 100 - round(100 * (total_records - passing_filter) / total_records, 2)
    print("Read {} total records, {} filtered records ({}%) in {} seconds".format(total_records, passing_filter, percent_passing, elapsed))
    return {'Total Raw Reads': total_records, 'Filtered Reads': passing_filter}

def unzip_files(sample=None, isCloud=False):
    # unzips files for a given sample
    # Files should be in the 'raw' directory, either in a folder starting
    # with the sample id, or a filename starting with the sample id
    raw_files_dir = Path(read_files_dir)
    if isCloud:
        raw_files_dir = os.path.join(Path(__file__).parent.parent.absolute(), 'tmp', 'raw')
    raw_path = Path(raw_files_dir)
    zip_files = list(chain(raw_path.glob(f'{sample}*.zip'), raw_path.glob(f'**/{sample}*/*.zip')))
    for zipped in zip_files:
        with zipfile.ZipFile(zipped.resolve(), 'r') as zip_ref:
            zip_ref.extractall(raw_files_dir)

    gzip_files = list(chain(raw_path.glob(f'{sample}*.gz'), raw_path.glob(f'**/{sample}*/*.gz')))
    for zipped in gzip_files:
        if Path(str(zipped.resolve())[:-3]).exists():
            continue
        with gzip.open(zipped.resolve()) as f_in:
            with open(str(zipped.resolve())[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                