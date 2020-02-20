import os
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
from pipeline.read_raw_files import filtergen

rawFilesDir = "../raw_files"
Qscore_threshold = 20

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
	print("Filtering records for {}...".format(code))
	filenames = [path.resolve() for path in Path(rawFilesDir).glob(code + '*.fastq')]
	total_records = 0

	concat_path = inter_path('{}_CONCAT.fastq'.format(code))
	with open(inter_path(concat_path), 'w') as outfile:
		for file in filenames:
			with open(file, 'r') as fastq_file:
				for line in fastq_file:
					outfile.write(line)
					if line.startswith('@'):
						total_records += 1

	filtered_records = 0
	threshold = Qscore_threshold

	filtered_path = inter_path('{}_{}_FILTERED.fastq'.format(code, threshold))
	with open(filtered_path, 'w+') as outfile:
		filtered_records += SeqIO.write(filtergen(concat_path, Qscore_threshold), outfile, "fastq")

	os.remove(concat_path)
	print("{} total records, {} filtered records for code {}".format(total_records, filtered_records, code))


