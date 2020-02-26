from pathlib import Path
import os
import csv
import fnmatch
import zipfile
import gzip
import shutil
from Bio import SeqIO
import itertools

from pipeline.utils import inter_path, update_log, output_path, setup_paths

from pipeline.read_raw_files import process_files, unzip_files
from pipeline.fingerprinting import fingerprinting
from pipeline.read_aligner import run_alignment
from pipeline.genome_histogram import make_plots

from parameters import raw_files_dir, Qscore_threshold, genome_path

def main():
	# First unzip any zip files (deletes the zips if "delete_intermediates" is true)
	unzip_files()

	# Examine all raw fastq files to see what codes are there for processing
	raw_files = Path(raw_files_dir).glob('*.fastq')
	codes = set(r.name.split('_')[0] for r in raw_files)
	for code in codes:
		print("PROCESSING CODE {}...".format(code))
		print('----------')
		print('----------')
		setup_paths(code, Qscore_threshold)

		log_info = {'Code': code, 'Qscore Threshold': str(Qscore_threshold)}
		run_prefix = "{}_Q{}".format(code, Qscore_threshold)
		
		# Start at the end to avoid repeating steps with saved results
		mapped_path = output_path("{}_unique_reads_aligned_histogram.csv".format(run_prefix))
		if not Path(mapped_path).exists():
			print("can't find {}".format(mapped_path))
			fingerprinted_path = inter_path("{}_FINGERPRINTED.fasta".format(run_prefix))
			if not Path(fingerprinted_path).exists():
				
				filtered_path = inter_path('{}_FILTERED.fastq'.format(run_prefix))
				if not Path(filtered_path).exists():
					filenames = [path.resolve() for path in Path(raw_files_dir).glob(code + '*.fastq')]
					process_results = process_files(code, filenames, filtered_path)
					log_info.update(process_results)
				
				fp_results = fingerprinting(filtered_path, fingerprinted_path)
				log_info.update(fp_results)

			alignment_results = run_alignment(fingerprinted_path, run_prefix)
			log_info.update(alignment_results)
		print("Final log info for {}: {}".format(code, log_info))
		update_log(log_info)

def check():
	hist_data_file = "{}_Q{}_unique_reads_aligned_histogram.csv".format('A4632', Qscore_threshold)
	hist_file_path = Path(output_path(hist_data_file))
	if (hist_file_path.exists()):
		# check if matches a pres-existing version
		with open(hist_file_path) as newF:
			readNew = csv.reader(newF)
			compareName = '../Example/A4632_Trans_Sites(genome mapping output).csv'
			compareName = Path(output_path("{}_Q{}_unique_reads_aligned_histogram_original.csv".format('A4632', Qscore_threshold)))
			with open(compareName, 'r') as orig:
				readOld = csv.reader(orig)
				readOld.__next__()
				readNew.__next__()
				newRow = readNew.__next__()
				for oldRow in readOld:
					if oldRow[1] != '0':
						if int(oldRow[0]) > int(newRow[0]):
							newRow = readNew.__next__()
						if newRow[0] != oldRow[0] or newRow[1] != oldRow[1]:
							print("Problem", oldRow, newRow)

def plot():
	print("plotting")
	hist_data_file = "{}_Q{}_unique_reads_aligned_histogram.csv".format('A4632', Qscore_threshold)
	hist_file_path = Path(output_path(hist_data_file))
	if (hist_file_path.exists()):
		make_plots(hist_file_path)

#main()
setup_paths('A4632', Qscore_threshold)
#check()
plot()