py_name = 'read_map_v2.py'

# v2 - does not analyze for spike-ins
# No longer allows for mismatch when looking at donor reads
import os
import datetime
import csv
import fnmatch
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import xlsxwriter
from pathlib import Path

# change directory based on where input files are
os.chdir('C:\\Users\\Leo Vo\\Desktop\\New_NGS_Dec_2019\\zipped_fastq\\fingerprinted')

map_length = 17 # length in bp of fingerprint
TSD = 5

##INGORES NON-UNIQUE MAPPING READS (Saves them to a separate text file)

# read in BL21DE3 RefSeq
for record in SeqIO.parse("genome.fasta", "fasta"):
    genome = record.seq.upper()  # remember to convert to upper case
    genome_rc = genome.reverse_complement()
len_genome = len(genome)
# add bases to the end from front to correct for sequences at the linear junctions
genome = genome + genome[:map_length-1]
genome_rc = genome_rc + genome_rc[:map_length-1]


excel_out_name = "read_map_log_BL21_ATW.xlsx"
excel_out_path = Path(excel_out_name)
if not excel_out_path.is_file():
    log = xlsxwriter.Workbook(excel_out_name)
    bold = log.add_format({'bold': True})
    percentage_format = log.add_format()
    percentage_format.set_num_format('0.00%')
    logsheet = log.add_worksheet("Fingerprint Log")
    logsheet.write(0,0, "File Codename", bold)
    logsheet.write(0, 1, "Total Fingerprints", bold)
    logsheet.write(0, 2, "Donor Reads", bold)
    logsheet.write(0, 3, "% Donor Reads", bold)
    logsheet.write(0, 4, "Genome-Mapping Reads", bold)
    logsheet.write(0, 5, "Unique Genome-Mapping Reads", bold)

# function to find all occurences of each read
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)  # use start += 1 to find overlapping matches

def read_map(filename, line, codename):
    donor_fp = "TGCTGAAACCTCAGGCA"
    total_fp = 0
    total_donor = 0
    total_genome = 0
    total_dup = 0
    tally = [0] * len_genome  # list for tallying locations of transposition
    dup = open('trans_sites\\dup_reads\\{}_non_unique_reads.txt'.format(codename), 'w', newline='')
    for record in SeqIO.parse(filename, "fasta"):
        total_fp += 1
        if record.seq == donor_fp:
            total_donor += 1
        else:
            # map to genome and tally
            locations = list(find_all(genome, record.seq))  # look in forward strand first
            locations_rc = list(find_all(genome_rc, record.seq))  # then look in reverse strand

            if len(locations) == 1 and len(locations_rc) == 0: # unique on fwd strand
                new_rec = SeqRecord(record.seq, id=record.id, name=record.name, description=record.description)
                yield new_rec
                total_genome += 1
                for i in locations:
                    site = i + map_length
                    if site > len_genome:  # correct for junction indexing
                        site = site - len_genome
                    tally[site - 1] += 1  # -1 to correct for indexing

            if len(locations_rc) == 1 and len(locations) == 0: # unique on rev strand
                new_rec = SeqRecord(record.seq, id=record.id, name=record.name, description=record.description)
                yield new_rec
                total_genome += 1
                for m in locations_rc:
                    site_rc = len_genome - m - map_length + TSD  # convert to fwd strand coordinates
                    if site_rc < 0:
                        site_rc = site_rc + len_genome
                    tally[site_rc - 1] += 1  # -1 to correct for indexing
            if len(locations) + len(locations_rc) > 1:
                total_dup += 1
                dup.write(str(record.seq)+'\n')
    with open('trans_sites\\{}_Trans_Sites.csv'.format(codename), 'w', newline='') as csv_out:
        writer = csv.writer(csv_out)
        writer.writerow(['Position', 'Coverage'])  # include header row
        for i in range(0, len(tally)):
            writer.writerow([i + 1, tally[i]])
    csv_out.close()
    dup.close()
    if os.stat('trans_sites\\dup_reads\\{}_non_unique_reads.txt'.format(codename)).st_size == 0:
        os.remove('trans_sites\\dup_reads\\{}_non_unique_reads.txt'.format(codename))
    logsheet.write(line, 1, total_fp)
    logsheet.write(line, 2, total_donor)
    logsheet.write(line, 3, total_donor / total_fp, percentage_format)
    logsheet.write(line, 4, total_genome + total_dup)
    logsheet.write(line, 5, total_genome)

print(datetime.datetime.now())

# read in the info csv file and begin processing files as a batch
with open('input_BL21_atw.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader, None)  # skips header row
    line = 1
    for row in reader:
        codename = row[0]
        print(codename)
        # search for file in cwd based on Sample ID ('code')
        filename = 'none'
        for i in os.listdir('.'):
            if fnmatch.fnmatch(i, "{}_FP.fasta".format(codename)):
                filename = i
                break
        if filename != 'none':
            logsheet.write(line, 0, codename)
            SeqIO.write(read_map(filename, line, codename), 'mapped_reads\\{}_unique.fasta'.format(codename), "fasta")
        else:
            print("WARNING - File Not Found For {}".format(codename))
        line += 1

log.close()  # close excel output
print(datetime.datetime.now())
