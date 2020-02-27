from pathlib import Path
import os
import csv
from tempfile import NamedTemporaryFile
import shutil
from datetime import datetime

intermediates_dir = ''
outputs_dir = ''
def setup_paths(code, Qscore):
	folder_name = '{}_Q{}'.format(code, Qscore)
	global intermediates_dir
	global outputs_dir
	intermediates_dir = os.path.join(Path(__file__).cwd(), 'intermediates', folder_name)
	outputs_dir = os.path.join(Path(__file__).cwd(), 'outputs', folder_name)
	os.makedirs(intermediates_dir, exist_ok=True)
	os.makedirs(os.path.join(outputs_dir, 'plots'), exist_ok=True)
	os.makedirs(outputs_dir, exist_ok=True)

def inter_path(filename):
	return os.path.join(intermediates_dir, filename)

def output_path(filename):
	return os.path.join(outputs_dir, filename)

log_fieldnames = [
	'Code', 
	'Genome Length',
	'Qscore Threshold',
	'Total Reads',
	'Filtered Reads',
	'Valid Fingerprint Reads',
	'Unique Genome-Mapping Reads',
	'Non-Unique Genome-Mapping Reads',
	'Donor Reads', 
	'Analysis Date'
]

def update_log(data):
	if not data['Code'] or not data['Qscore Threshold']:
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

			#for row in reader:
				#if row['Code'] == data['Code'] and row['Qscore Threshold'] == data['Qscore Threshold']:
				#	row.update(data)
				#	found = True
				#writer.writerow(row)
			if not found:
				writer.writerow(data)
		shutil.move(tempfile.name, log_path)

def get_log_entry(code, qscore):
	log_path = os.path.join(outputs_dir, '..', 'output_log.csv')
	if not Path(log_path).exists():
		print("Couldn't find the log file")
		return
	with open(log_path, 'r') as csvfile:
		reader = csv.DictReader(csvfile, fieldnames=log_fieldnames)
		for row in reader:
			if row['Code'] == code and row['Qscore Threshold'] == qscore:
				return row
