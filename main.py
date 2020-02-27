from pathlib import Path
import os
import csv
import fnmatch
import zipfile
import gzip
import shutil
from Bio import SeqIO
import itertools

from pipeline.utils import inter_path, update_log, output_path, setup_paths, get_info_for_sample

from pipeline.read_raw_files import process_files, unzip_files
from pipeline.fingerprinting import fingerprinting
from pipeline.read_aligner import run_alignment
from pipeline.plotting import make_genome_plots
from pipeline.trans_dist_plot import make_trans_dist_plot

from parameters import raw_files_dir, Qscore_threshold, genome_path, info_file

def main():
	if not Path(info_file).exists():
		print("Could not find the info file with sample information, see the parameters.py file")
		return

	# First unzip any zip files (deletes the zips if "delete_intermediates" is true)
	unzip_files()

	# Examine all raw fastq files to see what samples are there for processing
	raw_files = Path(raw_files_dir).glob('*.fastq')
	samples = set(r.name.split('_')[0] for r in raw_files)
	for sample in samples:
		print("PROCESSING SAMPLE {}...".format(sample))
		print('----------')
		print('----------')
		setup_paths(sample, Qscore_threshold)

		log_info = {'Sample': sample, 'Qscore Threshold': str(Qscore_threshold)}
		meta_info = get_info_for_sample(sample)

		run_prefix = "{}_Q{}".format(sample, Qscore_threshold)

		# Start at the end to avoid repeating steps with saved results
		histogram_path = output_path("{}_genome_read_locations.csv".format(run_prefix))
		unique_reads_path = output_path("{}_unique_reads.fasta".format(run_prefix))
		if not Path(histogram_path).exists() or not Path(unique_reads_path).exists():

			fp_path = inter_path("{}_FINGERPRINTED.fasta".format(run_prefix))
			if not Path(fp_path).exists():

				filtered_path = inter_path('{}_FILTERED.fastq'.format(run_prefix))
				if not Path(filtered_path).exists():
					filenames = [path.resolve() for path in Path(raw_files_dir).glob(sample + '*.fastq')]
					process_results = process_files(sample, filenames, filtered_path)
					log_info.update(process_results)
				
				fp_results = fingerprinting(filtered_path, fp_path)
				log_info.update(fp_results)

			alignment_results = run_alignment(fp_path, run_prefix, meta_info)
			log_info.update(alignment_results)

		update_log(log_info)
		run_information = make_genome_plots(histogram_path, meta_info)
		
		make_trans_dist_plot(unique_reads_path, run_information)
		
main()