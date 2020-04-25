from pathlib import Path
import os
import csv
import fnmatch
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
from pipeline.plasmid_plot import plot_plasmid

from parameters import working_dir, Qscore_threshold, info_file, delete_intermediates

def main():
	if not Path(info_file).exists():
		print("Could not find the info file with sample information, see the parameters.py file")
		return

	# First unzip any zip files (deletes the zips if "delete_intermediates" is true)
	unzip_files()

	samples_to_process = []
	with open(Path(info_file), 'r', encoding='utf-8-sig') as opened_info_file:
		reader = csv.DictReader(opened_info_file)
		samples_to_process = [row for row in reader]

	for sample_info in samples_to_process:
		sample = sample_info['Sample']
		meta_info = get_info_for_sample(sample)
		if not meta_info:
			print("Sample {} not found in the info file, skipping processing".format(sample))
			continue

		print('----------')
		print('----------')
		print("PROCESSING SAMPLE {}...".format(sample))
		print('----------')
		print('----------')
		setup_paths(sample)

		log_info = {'Sample': sample, 'Qscore Threshold': str(Qscore_threshold)}

		# Start at the end to avoid repeating steps with saved results
		histogram_path = output_path(os.path.join('samples', "{}_genome_read_locations.csv".format(sample)))
		plasmid_histogram_path = output_path(os.path.join('samples', "{}_plasmid_read_locations.csv".format(sample)))
		
		filtered_path = inter_path('{}_FILTERED.fastq'.format(sample))
		fp_path = inter_path("{}_FINGERPRINTED.fasta".format(sample))

		# step 1: process raw files, concatenate
		raw_files_dir = os.path.join(Path(working_dir), 'raw')
		filtered_path = inter_path('{}_FILTERED.fastq'.format(sample))
		filenames = [path.resolve() for path in Path(raw_files_dir).glob(sample + '*.fastq')]
		if len(filenames) < 1:
			print("COULD NOT FIND ANY FASTA FILES FOR THE SAMPLE")
			continue
		process_results = process_files(sample, filenames, filtered_path)
		log_info.update(process_results)

		# step 2: fingerprint the reads
		fp_path = inter_path("{}_FINGERPRINTED.fasta".format(sample))
		fp_results = fingerprinting(filtered_path, fp_path, meta_info)
		log_info.update(fp_results)

		# step 3: align the reads in the genome
		alignment_results = run_alignment(fp_path, meta_info)
		log_info.update(alignment_results)

		update_log(log_info)
		if delete_intermediates:
			shutil.rmtree(inter_path(''))

		run_information = make_genome_plots(histogram_path, meta_info)

		if len(meta_info['Plasmid fasta file']) > 1:
			run_information = plot_plasmid(plasmid_histogram_path, meta_info)
		
		make_trans_dist_plot(histogram_path, run_information)
	
	if delete_intermediates:
		shutil.rmtree(Path(inter_path('')).parent.absolute())

main()