from pathlib import Path
import os
import csv
from tempfile import NamedTemporaryFile
import shutil
from datetime import datetime

from parameters import info_file, working_dir

intermediates_dir = ''
outputs_dir = ''
def setup_paths(code, isCloud):
	global intermediates_dir
	global outputs_dir
	repo_dir = Path(__file__).parent.parent.absolute()
	working = working_dir
	if isCloud:
		working = os.path.join(repo_dir, 'tmp')
		os.makedirs(os.path.join(working, 'raw'), exist_ok=True)
	intermediates_dir = os.path.join(Path(working), 'intermediates', code)
	outputs_dir = os.path.join(repo_dir, 'outputs')
	os.makedirs(intermediates_dir, exist_ok=True)
	os.makedirs(os.path.join(outputs_dir, 'samples'), exist_ok=True)
	os.makedirs(os.path.join(outputs_dir, 'plots'), exist_ok=True)
	os.makedirs(outputs_dir, exist_ok=True)

def inter_path(filename):
	return os.path.join(intermediates_dir, filename)

def output_path(filename):
	return os.path.join(outputs_dir, filename)

log_fieldnames = [
	'Sample', 
	'Qscore Threshold',
	'Total Raw Reads',
	'Filtered Reads',
	'Valid Fingerprint Reads',
	'Unique Genome-Mapping Reads',
	'Total Genome-Mapping Reads',
	'Undigested Donor Reads', 
	'Spike-in Reads',
	'CRISPR Array Self-Targeting Reads',
	'Analysis Date',
	'Unique Plasmid-Mapping Reads',
	'Total Plasmid-Mapping Reads'
]

def update_log(data):
	if not data['Sample'] or not data['Qscore Threshold']:
		print("Cannot make a log entry without a code and Qscore threshold")
		return

	data['Analysis Date'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

	log_path = os.path.join(outputs_dir, 'output_log.csv')
	if not Path(log_path).exists():
		with open(log_path, 'w') as csvfile:
			writer = csv.DictWriter(csvfile, fieldnames=log_fieldnames, extrasaction='ignore')
			writer.writeheader()
			writer.writerow(data)

	else:
		tempfile = NamedTemporaryFile(mode='w', delete=False)
		with open(log_path, 'r') as csvfile, tempfile:
			reader = csv.DictReader(csvfile, fieldnames=log_fieldnames)
			writer = csv.DictWriter(tempfile, fieldnames=log_fieldnames, extrasaction='ignore')
			found = False

			for row in reader:
				# Use this if you want to update existing rows instead of add new ones
				if row['Sample'] == data['Sample'] and row['Qscore Threshold'] == data['Qscore Threshold']:
					row.update(data)
					found = True
					
				writer.writerow(row)
			if not found:
				writer.writerow(data)
		shutil.move(tempfile.name, log_path)

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
