from pathlib import Path
import os
import csv
from tempfile import NamedTemporaryFile
import shutil
from datetime import datetime

from parameters import info_file, read_files_dir

intermediates_dir = ''
outputs_dir = ''
def setup_paths(code, isCloud):
	global intermediates_dir
	global outputs_dir
	repo_dir = Path(__file__).parent.parent.absolute()
	if isCloud:
		working = os.path.join(repo_dir, 'tmp')
		os.makedirs(os.path.join(working, 'raw'), exist_ok=True)
	intermediates_dir = os.path.join(repo_dir, 'intermediates', code)
	outputs_dir = os.path.join(repo_dir, 'outputs')
	os.makedirs(intermediates_dir, exist_ok=True)
	os.makedirs(os.path.join(outputs_dir, 'samples'), exist_ok=True)
	os.makedirs(os.path.join(outputs_dir, 'plots'), exist_ok=True)
	os.makedirs(outputs_dir, exist_ok=True)

def inter_path(filename):
	if filename:
		return os.path.join(intermediates_dir, filename)
	return intermediates_dir

def output_path(filename):
	return os.path.join(outputs_dir, filename)

log_fieldnames = [
	'Sample', 
	'Experiment date',
	'Run date',
	'Qscore Threshold',
	'Total Raw Reads',
	'Filtered Reads',
	'Transposon end-containing reads',
	'Unique Target Mapping Reads',
	'Total Target Mapping Reads',
	'Contaminating donor reads', 
	'Spike-in Reads',
	'CRISPR Array Self-Targeting Reads',
	'Analysis Date',
	'Unique Second Target Mapping Reads',
	'Total Second Target Mapping Reads'
]

def update_log(all_log_data):
	analysis_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

	log_path = os.path.join(outputs_dir, f'{datetime.now().strftime("%Y%m%d-%H:%M")}_output_log.csv')
	fieldnames = log_fieldnames + [
		'% Transposon end-containing',
		'% Unique target',
		'% Unique second target',
		'% Donor',
		'% Spike-in',
		'% CRISPR Array',
		'% Other'
	]
	with open(log_path, 'w') as csvfile:
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction='ignore')
		writer.writeheader()
		for data in all_log_data:
			data['Analysis Date'] = analysis_date
			data['% Transposon end-containing'] = int(data['Transposon end-containing reads']) / int(data['Filtered Reads'])
			data['% Unique target'] = int(data['Unique Target Mapping Reads']) / int(data['Filtered Reads'])
			data['% Unique second target'] = int(data.get('Unique Second Target Mapping Reads', 0)) / int(data['Filtered Reads'])
			data['% Donor'] = int(data['Contaminating donor reads']) / int(data['Filtered Reads'])
			data['% Spike-in'] = int(data['Spike-in Reads']) / int(data['Filtered Reads'])
			data['% CRISPR Array'] = int(data['CRISPR Array Self-Targeting Reads']) / int(data['Filtered Reads'])
			data['% Other'] = 1 - sum([
				data['% Unique target'], 
				data['% Unique second target'],
				data['% Donor'],
				data['% Spike-in'],
				data['% CRISPR Array']
			])
			writer.writerow(data)
	return log_path

def get_log_entry(sample, qscore):
	log_path = os.path.join(outputs_dir, '..', 'output_log.csv')
	return get_row_from_csv(sample, log_path, qscore)

def get_info_for_sample(sample):
	return get_row_from_csv(sample, info_file)

def get_row_from_csv(sample, filepath, qscore=None):
	if not Path(filepath).exists():
		print(f"Couldn't find the file {filepath}")
		return
	result = None
	try:
		with open(Path(filepath), 'r', encoding='utf-8-sig') as csvfile:
			reader = csv.DictReader(csvfile)
			for row in reader:
				if row['Sample'] == sample:
					if not qscore or row['Qscore Threshold'] == qscore:
						result = row
						break
	except:
		with open(Path(filepath), 'r', encoding='ISO-8859-1') as csvfile:
			reader = csv.DictReader(csvfile)
			for row in reader:
				if row['Sample'] == sample:
					if not qscore or row['Qscore Threshold'] == qscore:
						result = row
						break
	return result
