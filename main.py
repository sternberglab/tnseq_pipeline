from pathlib import Path
import os
import csv
import fnmatch
import zipfile
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import itertools

from pipeline.utils import inter_path

from pipeline.read_raw_files import process_files
from pipeline.fingerprinting import fingerprinting

from parameters import rawFilesDir, Qscore_threshold

# First unzip any zip files
zip_files = Path(rawFilesDir).glob('*.zip')
for zipped in zip_files:
	with zipfile.ZipFile(zipped.resolve(), 'r') as zip_ref:
		zip_ref.extractall(rawFilesDir)

# Get all raw fastq files
# Do the initial filter to an intermediate file
raw_files = Path(rawFilesDir).glob('*.fastq')
codes = set(r.name.split('_')[0] for r in raw_files)
for code in codes:
	print("Processing code {}".format(code))
	
	filtered_path = inter_path('{}_Q{}_FILTERED.fastq'.format(code, Qscore_threshold))
	if not Path(filtered_path).exists():
		print("Filtering records for {}...".format(code))
		filenames = [path.resolve() for path in Path(rawFilesDir).glob(code + '*.fastq')]
		total_records, filtered_records = process_files(filenames, filtered_path)
		print("{} total records, {} filtered records for code {} - ({}%)".format(total_records, filtered_records, code, filtered_records/total_records))

	fingerprinted_path = inter_path("{}_Q{}_FINGERPRINTED.fastq".format(code, Qscore_threshold))
	if not Path(fingerprinted_path).exists():
		print("Fingerprinting records for {}...".format(code))
		total_reads, fp_reads = fingerprinting(filtered_path, fingerprinted_path)
		print("Of {} reads, {} ({}%) had a valid fingerprint".format(total_reads, fp_reads, fp_reads/total_reads))



