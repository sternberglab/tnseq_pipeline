from pathlib import Path
import os
import copy
import csv
import json
import fnmatch
import gzip
import shutil
from Bio import SeqIO
import itertools
import botocore
import boto3

from pipeline.utils import inter_path, update_log, output_path, setup_paths, get_info_for_sample

from pipeline.read_raw_files import process_files, unzip_files
from pipeline.fingerprinting import fingerprinting
from pipeline.read_aligner import run_alignment
from pipeline.plotting import make_genome_plots
from pipeline.trans_dist_plot import make_trans_dist_plot
from pipeline.plasmid_plot import plot_plasmid

from parameters import working_dir, Qscore_threshold, info_file, delete_intermediates

def get_samples_to_process(isCloud):
	if not isCloud:
		# if running locally use this path
		if not Path(info_file).exists():
			print("Could not find the info file with sample information, see the parameters.py file")
			return
		samples_to_process = []
		try:
			with open(Path(info_file), 'r', encoding='utf-8-sig') as opened_info_file:
				reader = csv.DictReader(opened_info_file)
				samples_to_process = [row for row in reader]
		except: 
			with open(Path(info_file), 'r', encoding='ISO-8859-1') as opened_info_file:
				reader = csv.DictReader(opened_info_file)
				samples_to_process = [row for row in reader]
		return samples_to_process

	# if on aws load from the sqs json and download the data files from s3
	with open('./sqs_message.json', 'r') as info_json:
		info = json.load(info_json)
	setup_paths(info['Sample'], isCloud)
	# if in cloud download the s3 files to the container
	if isCloud:
		s3 = boto3.client('s3')
		for item in s3.list_objects(Bucket='sternberg-sequencing-data', Prefix=f"ngs_data/basespace/{info['analysisId']}/{info['Sample']}")['Contents']:
			key = item['Key']
			filename = key.split('/')[-1]
			s3.download_file('sternberg-sequencing-data', key, f'./tmp/raw/{filename}')
		if info['Genome']:
			s3.download_file('sternberg-sequencing-data', f"bioinformatic_resources/genomes/{info['Genome']}", f"./tmp/{info['Genome']}")
		if info['Plasmid']:
			s3.download_file('sternberg-sequencing-data', f"bioinformatic_resources/plasmids/{info['Plasmid']}", f"./tmp/{info['Plasmid']}")

	return [info]

def main(isCloud=False):
	samples_to_process = get_samples_to_process(isCloud)
	for sample_info in samples_to_process:
		sample = sample_info['Sample']
		# unzip files for the sample (deletes the zips if "delete_intermediates" is true)
		unzip_files(sample, isCloud)
		setup_paths(sample, isCloud)
		meta_info = sample_info if isCloud else get_info_for_sample(sample)
		original_input = copy.deepcopy(meta_info)
		if isCloud:
			meta_info['Genome'] = os.path.join(Path(__file__).parent.absolute(), 'tmp', meta_info['Genome'])
			if meta_info['Plasmid']:
				meta_info['Plasmid'] = os.path.join(Path(__file__).parent.absolute(), 'tmp', meta_info['Plasmid'])

		print('----------')
		print('----------')
		print("PROCESSING SAMPLE {}...".format(sample))
		print('----------')
		print('----------')


		log_info = {'Sample': sample, 'Qscore Threshold': str(Qscore_threshold), 'Input Parameters': original_input}

		
		# Start at the end to avoid repeating steps with saved results
		histogram_path = output_path(os.path.join('samples', "{}_genome_read_locations.csv".format(sample)))
		plasmid_histogram_path = output_path(os.path.join('samples', "{}_plasmid_read_locations.csv".format(sample)))
		
		filtered_path = inter_path('{}_FILTERED.fastq'.format(sample))
		fp_path = inter_path("{}_FINGERPRINTED.fasta".format(sample))
		# step 1: process raw files, concatenate
		raw_files_dir = os.path.join(Path(working_dir), 'raw') if not isCloud else './tmp/raw'
		filtered_path = inter_path('{}_FILTERED.fastq'.format(sample))
		files = list(itertools.chain(Path(raw_files_dir).glob(f"{sample}*.fastq"), Path(raw_files_dir).glob(f"{sample}*/*.fastq")))
		filenames = [path.resolve() for path in files]
		if len(filenames) < 1:
			print("COULD NOT FIND ANY FASTA FILES FOR THE SAMPLE")
			continue
		process_results = process_files(sample, filenames, filtered_path, meta_info)
		log_info.update(process_results)
		
		# step 2: fingerprint the reads
		fp_path = inter_path("{}_FINGERPRINTED.fasta".format(sample))
		fp_results = fingerprinting(filtered_path, fp_path, meta_info)
		log_info.update(fp_results)

		# step 3: align the reads in the genome
		alignment_results = run_alignment(fp_path, meta_info)
		log_info.update(alignment_results)

		update_log(log_info)

		run_information = make_genome_plots(histogram_path, meta_info)

		if len(meta_info['Plasmid']) > 1:
			run_information = plot_plasmid(plasmid_histogram_path, meta_info)
		
		make_trans_dist_plot(histogram_path, run_information)
		
		if delete_intermediates:
			# delete intermediate files
			shutil.rmtree(Path(inter_path('')).absolute())
			# delete unzipped files if zips are present
			for file in files:
				if Path(f'{file}.gz').exists():
					os.remove(file)

		
	return log_info

if __name__== "__main__":
	main()