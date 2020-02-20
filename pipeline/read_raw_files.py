py_name = 'read_raw_files.py'


###############
# Script Function: 3-in-1 script
# 1/ Extracts .gz files from individual subfolders to return fastq files (usually each fastq correspond to 1 lane of 4)
# 2/ Concatenate all separate fastq files for lanes (usually 4) of one sample into 1 main fastq per sample
# 3/ Reads with more than half bases with low Q score are removed.
# (v2 - no longer switches low quality bases to N to recover more potential reads)
# 'Low' Qscore is determined by a threshold variable, and function returns a new fastq files with edited filtered reads.
# Keeps track of filtering stats with an excel spreadsheet output

# Script Input(s):
# 1/ One batch info file (CSV) containing numbers 1 to number of files in column 1,
# and Sample IDs (or 'codes') e.g. 'A4715' in column 2. Do not remove header row
# 2/ If running functions separately, note that concat and filtergen fucntions take .fastq files, which are usually
# outputs from previous functions in the script. Functions will output warnings accordingly if input files are not found
# 3/ NEED 1 SUBFOLDER CALLED 'IndvFastQs'

# Script Output(s):
# See individual functions. If run all together, returns one filtered fastq file for each SampleID/code/codename,
# and one main excel output for stats keeping purposes.
###############

import os
import fileinput
import csv
import fnmatch
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import xlsxwriter
import pathlib
from pathlib import Path
from .utils import inter_path


def filtergen(file, threshold):  # generator function that returns edited reads that pass filter, to write new fastq file
    for record in SeqIO.parse(file, "fastq"):
        # Convert base qualities to Boolean based on Qscore threshold value. Only use reads with >=50% non-N:
        recordqual = [x > threshold for x in record.letter_annotations['phred_quality']]  # list of True, False etc
        if float(sum(recordqual)) / float(len(recordqual)) >= .5:  # note that True = 1, False = 0 for summing
            # create new SeqRecord with edited read sequence
            newrec = SeqRecord(record.seq, id=record.id, name=record.name,
                               description=record.description, letter_annotations=record.letter_annotations)
            yield newrec