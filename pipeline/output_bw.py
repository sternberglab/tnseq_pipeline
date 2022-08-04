# This file creates IGV-readable xml files
# from samples and genomes referenced in the 
# input info_file (as set in parameters)
# 
# It requires the outputs of the pipeline being in their default
# locations in the outputs folder of the Illumina-pipeline directory


####### Actual code here ######
###############################
###############################
import os
from pathlib import Path
import subprocess
from parameters import info_file
from pipeline.utils import inter_path
import csv
from Bio import SeqIO

directory = Path(__file__).parent.parent.absolute()
igv_dir = os.path.join(directory, 'intermediates', 'igv')

def get_input_info():
	# Get the sample's input information
	samples_inputs = None
	info_file_path = Path(info_file)
	try:
		with open(info_file_path, 'r', encoding='utf-8-sig') as opened_info_file:
			reader = csv.DictReader(opened_info_file)
			all_inputs_info = [row for row in reader]
	except: 
		with open(info_file_path, 'r', encoding='ISO-8859-1') as opened_info_file:
			reader = csv.DictReader(opened_info_file)
			all_inputs_info = [row for row in reader]
	return {s["Sample"]: s for s in all_inputs_info}

def get_output_info(output_path):
	# Get the samples output information
	samples_outputs = None
	if not Path(output_path).exists():
		raise Exception(f"Could not find {output_path} to create bigwig files")
	with open(output_path, 'r', encoding='utf-8-sig') as opened_output_file:
		reader = csv.DictReader(opened_output_file)
		samples_outputs = {row["Sample"]: row for row in reader}
	return samples_outputs

def get_chromosome(filepath):
	path = Path(filepath)
	# gets the chromosome name - whatever is after ">" in the "Target fasta file"
	try:
		chromosome = SeqIO.read(path, 'fasta')
	except:
		try:
			chromosome = SeqIO.read(os.path.join(directory, path), 'fasta')
		except:
			raise Exception(f"Could not find the target fasta file at {filepath}")
	return chromosome

def create_bed_file(input_info, output_info, chromosome):
	sample_id = input_info["Sample"]

	# Get the chromosome name - whatever is after ">" in the "Target fasta file"
	chromosome_id = chromosome.id

	# Find the target_read_locations.csv file in the outputs folder
	output_path_date = output_info.get("Experiment date", None)
	if not output_path_date:
		output_path_date = output_info.get("Run date", None)
	target_reads_csv_path = os.path.join(directory, 'outputs', 'samples', f"{output_path_date}_{sample_id}_target_read_locations.csv")
	with open(target_reads_csv_path, 'r', encoding='utf-8-sig') as opened_reads_file:
		reader = csv.DictReader(opened_reads_file)
		# Create the new csv file for writing
		new_bed_filepath = os.path.join(igv_dir, sample_id, f"1_{sample_id}_all.bed")
		with open(new_bed_filepath, 'w') as sample_bed_out:
			writer = csv.writer(sample_bed_out, delimiter='\t')
			for row in reader:
				start = row['position']
				end = str(int(start)+1)
				# Row values should be:
				# 1. Chromosome id - whatever is after ">" in the "Target fasta file"
				# 2. Start - reference position of the read
				# 3. End - start + 1
				# 4. Total - the total reads from the output reads location csv
				# 5. RL - from the output reads location csv
				# 6. LR - from the output reads location csv
				writer.writerow([chromosome_id, start, end, row['reads'], row['RL'], row['LR']])
	return new_bed_filepath

def process_bed_to_bigwig(sample_id, bed_path,  output_info, chromosome):
	sample_dir = os.path.join(igv_dir, sample_id)

	# Make 3 bedgraph files from the csv, for total, lr, and rl
	total_path = os.path.join(sample_dir, f"2_parsed_{sample_id}_total.bedgraph")
	cmd_total = f"""cut -f1-4 {bed_path} > {total_path};"""
	subprocess.run(cmd_total, shell=True)

	rl_path = os.path.join(sample_dir, f"2_parsed_{sample_id}_RL.bedgraph")
	cmd_rl = f"""cut -f1-3,5 {bed_path} > {rl_path};"""
	subprocess.run(cmd_rl, shell=True)
	
	lr_path = os.path.join(sample_dir, f"2_parsed_{sample_id}_LR.bedgraph")
	cmd_lr = f"""cut -f1-3,6 {bed_path} > {lr_path};"""
	subprocess.run(cmd_lr, shell=True)

	# Create chrom.sizes file
	chrom_sizes_path = os.path.join(sample_dir, 'chrom.sizes')
	with open(chrom_sizes_path, 'w') as sizes_file:
		sizes_file.write(f"{chromosome.id}\t{len(chromosome.seq)}")
	
	# Run bedGraphToBigWig on each
	outputs_igv_dir = os.path.join(directory, 'outputs', 'igv', sample_id)
	os.makedirs(outputs_igv_dir, exist_ok=True)
	for ct_type in ['total', 'RL', 'LR']:
		input_filepath = os.path.join(sample_dir, f"2_parsed_{sample_id}_{ct_type}.bedgraph")
		
		# Output to the outputs directory, not intermediates
		# Make the output directory for igv if not exists
		date_prefix = output_info['Experiment date']
		if not date_prefix:
			date_prefix = output_info["Run date"]
		output_filepath = os.path.join(outputs_igv_dir, f'{date_prefix}_{sample_id}_{ct_type}.bw')
		to_bigwig_cmd = f"bedGraphToBigWig {input_filepath} {chrom_sizes_path} {output_filepath}"
		subprocess.run(to_bigwig_cmd, shell=True)

def create_bw_outputs(output_path):
	# Get sample info from the inputs and outputs
	all_input_info = get_input_info()
	all_output_info = get_output_info(output_path)
	
	# For each sample, create a csv that can be converted to bed format
	# and place in the "1_all_csv" folder
	for sample_id in all_input_info.keys():
		os.makedirs(os.path.join(igv_dir, sample_id), exist_ok=True)
		input_info = all_input_info[sample_id]
		output_info = all_output_info.get(sample_id, None)
		if not output_info or int(output_info['Unique Target Mapping Reads']) == 0:
			continue

		chromosome = get_chromosome(input_info["Target fasta file"])
		bed_path = create_bed_file(input_info, output_info, chromosome)
		process_bed_to_bigwig(sample_id, bed_path, output_info, chromosome)
