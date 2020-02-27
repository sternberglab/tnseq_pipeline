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
from pipeline.plotting import make_genome_plots
from pipeline.trans_dist_plot import make_trans_dist_plot

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
		histogram_path = output_path("{}_genome_read_locations.csv".format(run_prefix))
		fp_path = output_path("{}_FINGERPRINTED.fasta".format(run_prefix))
			
		if not Path(histogram_path).exists() or not Path(fp_path).exists():
			print("can't find either histogram csv or fingerprinted reads for {}".format(run_prefix))
			
			if not Path(fp_path).exists():
				filtered_path = inter_path('{}_FILTERED.fastq'.format(run_prefix))
				if not Path(filtered_path).exists():
					filenames = [path.resolve() for path in Path(raw_files_dir).glob(code + '*.fastq')]
					process_results = process_files(code, filenames, filtered_path)
					log_info.update(process_results)
				
				fp_results = fingerprinting(filtered_path, fp_path)
				log_info.update(fp_results)

			alignment_results = run_alignment(fp_path, run_prefix)
			log_info.update(alignment_results)

		update_log(log_info)

		run_information = make_genome_plots(histogram_path)
		
		make_trans_dist_plot(fp_path, run_information)
		

main()