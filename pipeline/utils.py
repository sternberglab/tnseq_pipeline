from pathlib import Path
import os
import csv
from tempfile import NamedTemporaryFile
import shutil
from datetime import datetime

from parameters import info_file, working_dir

intermediates_dir = ''
outputs_dir = ''
def setup_paths(code, Qscore):  # Set up path and folders for intermediate and output files
	folder_name = '{}_Q{}'.format(code, Qscore)
	global intermediates_dir
	global outputs_dir
	intermediates_dir = os.path.join(Path(working_dir), 'intermediates', folder_name)
	outputs_dir = os.path.join(Path(working_dir), 'outputs', folder_name)
	os.makedirs(intermediates_dir, exist_ok=True)
	os.makedirs(os.path.join(outputs_dir, 'plots'), exist_ok=True)
	os.makedirs(outputs_dir, exist_ok=True)

def inter_path(filename):  # Make a path in the intermediate files dir for a specified filename
	return os.path.join(intermediates_dir, filename)

def output_path(filename):  # Make a path in the output files dir for a specified filename
	return os.path.join(outputs_dir, filename)

# Types of analysis data stored into the log .csv
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
	'Analysis Date'
]

def update_log(data):  # writes analysis statistics into csv log
	if not data['Sample'] or not data['Qscore Threshold']:
		print("Cannot make a log entry without a code and Qscore threshold")
		return

	data['Analysis Date'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

	log_path = os.path.join(outputs_dir, '..', 'output_log.csv')
	if not Path(log_path).exists():
		with open(log_path, 'w') as csvfile:
			writer = csv.DictWriter(csvfile, fieldnames=log_fieldnames)
			writer.writeheader()
			writer.writerow(data)

	else:
		tempfile = NamedTemporaryFile(mode='w', delete=False)
		with open(log_path, 'r') as csvfile, tempfile:
			reader = csv.DictReader(csvfile, fieldnames=log_fieldnames)
			writer = csv.DictWriter(tempfile, fieldnames=log_fieldnames)
			found = False

			for row in reader:
				# Use this if you want to update existing rows instead of add new ones

				#if row['Sample'] == data['Sample'] and row['Qscore Threshold'] == data['Qscore Threshold']:
				#	row.update(data)
				#	found = True
				writer.writerow(row)
			if not found:
				writer.writerow(data)
		shutil.move(tempfile.name, log_path)

# def get_log_entry(sample, qscore):
# 	log_path = os.path.join(outputs_dir, '..', 'output_log.csv')
# 	if not Path(log_path).exists():
# 		print("Couldn't find the log file")
# 		return
# 	with open(log_path, 'r') as csvfile:
# 		reader = csv.DictReader(csvfile, fieldnames=log_fieldnames)
# 		for row in reader:
# 			if row['Sample'] == sample and row['Qscore Threshold'] == qscore:
# 				return row

def get_info_for_sample(sample):  # returns a dict of sample information for each sample from input csv
	with open(Path(info_file), 'r', encoding='utf-8-sig') as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			if row['Sample'] == sample:
				return row
	return
