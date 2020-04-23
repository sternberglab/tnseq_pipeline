from pathlib import Path
import os
import csv
import fnmatch
import gzip
import shutil
from Bio import SeqIO
import itertools

from .utils import inter_path, update_log, output_path, setup_paths, get_info_for_sample

from .read_raw_files import process_files, unzip_files
from .fingerprinting import fingerprinting
from .read_aligner import run_alignment
from .plotting import make_genome_plots
from .trans_dist_plot import make_trans_dist_plot
from .plasmid_plot import plot_plasmid

from parameters import working_dir, Qscore_threshold, info_file, delete_intermediates

def main(info_file, genome, system):
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
		meta_info = = sample_info

		print('----------')
		print('----------')
		print("PROCESSING SAMPLE {}...".format(sample))
		print('----------')
		print('----------')
		setup_paths(sample, Qscore_threshold)

		log_info = {'Sample': sample, 'Qscore Threshold': str(Qscore_threshold)}

		meta_info['run_prefix'] = sample

		# Start at the end to avoid repeating steps with saved results
		histogram_path = output_path("genome_read_locations.csv")
		plasmid_histogram_path = output_path("plasmid_read_locations.csv")
		
		unique_reads_path = output_path("genome_unique_reads.fasta")
		filtered_path = inter_path('{}_FILTERED.fastq'.format(run_prefix))
		fp_path = inter_path("{}_FINGERPRINTED.fasta".format(run_prefix))
		if not Path(histogram_path).exists() or not Path(plasmid_histogram_path).exists() or not Path(unique_reads_path).exists():
			fp_path = inter_path("{}_FINGERPRINTED.fasta".format(run_prefix))
			if not Path(fp_path).exists():
				filtered_path = inter_path('{}_FILTERED.fastq'.format(run_prefix))
				if not Path(filtered_path).exists():
					raw_files_dir = os.path.join(Path(working_dir), 'raw')
					filenames = [path.resolve() for path in Path(raw_files_dir).glob(sample + '*.fastq')]
					process_results = process_files(sample, filenames, filtered_path)
					log_info.update(process_results)
				fp_results = fingerprinting(filtered_path, fp_path)
				log_info.update(fp_results)
			alignment_results = run_alignment(fp_path, meta_info)
			update_log(log_info)
			if delete_intermediates:
				shutil.rmtree(os.path.join(Path(working_dir), 'intermediates', run_prefix))

		run_information = make_genome_plots(histogram_path, meta_info)

		if len(meta_info['Plasmid fasta file']) > 1:
			run_information = plot_plasmid(plasmid_histogram_path, meta_info)
		
		make_trans_dist_plot(unique_reads_path, run_information)