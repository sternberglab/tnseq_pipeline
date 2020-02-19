py_name = 'read_fingerprinting_v1.py'


import os
import csv
import fnmatch
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import xlsxwriter
from pathlib import Path

inputDir = '/Users/Chris_Acree/Documents/lab/Example/Pre-concatnation-fastqs'
#'C:\\Users\\Leo Vo\\Desktop\\New_NGS_Dec_2019\\zipped_fastq'

# change directory based on where input files are
os.chdir(inputDir)

excel_out_name = "fingerprint_log_2_122219.xlsx"
excel_out_path = Path(excel_out_name)
if not excel_out_path.is_file():
    log = xlsxwriter.Workbook(excel_out_name)
    bold = log.add_format({'bold': True})
    percentage_format = log.add_format()
    percentage_format.set_num_format('0.00%')
    logsheet = log.add_worksheet("Fingerprint Log")
    logsheet.write(0,0, "File Codename", bold)
    logsheet.write(0, 1, "Total Reads", bold)
    logsheet.write(0, 2, "Reads with valid Fingerprint", bold)
    logsheet.write(0, 3, "Percentage", bold)

def hamming_dist(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def fpgen(filename, flank, fp_length):  # generator function that returns
    for record in SeqIO.parse(filename, "fastq"):
        i = 0
        while 0 <= i <= (len(record.seq) - len(flank)):
            query = record.seq[i:i+len(flank)]
            if hamming_dist(query, flank) < 2:
                if i >= fp_length:
                    fingerprint = record.seq[i - fp_length:i]
                    newrec = SeqRecord(fingerprint, id=record.id, name=record.name)
                    yield newrec
                i += 2000
            else:
                i += 1

def fingerprinting(filename, line):  # main command function for filtering
    flank = "TGTTGGAACAACCAT"
    fp_length = 17
    total = 0
    for record in SeqIO.parse(filename, "fastq"):
        total += 1
    # generate input records to write new fingerprinted fastq
    fasta_filename = './fingerprinted/{}_FP.fasta'.format(codename)
    with open(fasta_filename, 'w+') as fastafile:
        SeqIO.write(fpgen(filename, flank, fp_length), fastafile, "fasta")
    fp_total = 0
    for record in SeqIO.parse(fasta_filename, "fasta"):
        fp_total += 1
    logsheet.write(line, 1, total)
    logsheet.write(line, 2, fp_total)
    logsheet.write(line, 3, fp_total / total, percentage_format)
    return

# read in the info csv file and begin processing files as a batch
with open('../../Example/input.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader, None)  # skips header row
    line = 1
    for row in reader:
        codename = row[0]
        # search for file in cwd based on Sample ID ('code')
        filename = 'none'
        for i in os.listdir('.'):
            if fnmatch.fnmatch(i, "{}_FINAL.fastq".format(codename)):
                filename = i
                break
        if filename != 'none':
            logsheet.write(line, 0, codename)
            fingerprinting(filename, line)
        else:
            print("WARNING - File Not Found For {}".format(codename))
        line += 1

log.close()  # close excel output

