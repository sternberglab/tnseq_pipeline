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
from pipeline.read_aligner import run_alignment

from parameters import rawFilesDir, Qscore_threshold, genome_path


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
	print("PROCESSING CODE {}...".format(code))
	print('----------')
	print('----------')
	log_info = {}
	run_prefix = "{}_Q{}".format(code, Qscore_threshold)
	
	# Start at the end to avoid repeating steps with saved results
	mapped_path = inter_path("{}_mapped.csv")
	if not Path(mapped_path).exists():
		
		fingerprinted_path = inter_path("{}_FINGERPRINTED.fasta".format(run_prefix))
		if not Path(fingerprinted_path).exists():
			
			filtered_path = inter_path('{}_FILTERED.fastq'.format(run_prefix))
			if not Path(filtered_path).exists():
				filenames = [path.resolve() for path in Path(rawFilesDir).glob(code + '*.fastq')]
				process_results = process_files(code, filenames, filtered_path)
				log_info.update(process_results)
			
			fp_results = fingerprinting(filtered_path, fingerprinted_path)
			log_info.update(fp_results)

		alignment_results = run_alignment(fingerprinted_path, run_prefix)
		log_info.update(alignment_results)
	print("Final log info for {}: {}".format(code, log_info))

