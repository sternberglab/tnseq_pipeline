This repository contains the code for interpreting and analyzing Illumina sequencing data for mapping transposon insertion events. 


# Getting started:

To run these scripts, you need to have Anaconda installed. Install this project by navigating to this directory and running `conda env create -n your_env_name --file environment.yml`, replacing `your_env_name` with your chosen environment name. 

# Running this code:

- Activate the conda environment: `conda activate your_env_name`. 
- Run the pipeline: `python main.py`

All necessary parameters should be set in `parameters.py` and in your input `.csv` file. 

# System parameters:

General settings are in the `parameters.py` file. See the comments on relevant parameters there, but at a minimum, you must set the following:

- `info_file`: A path to the csv file with information about the samples to run. See an example in './inputs/experiment_info.csv', and more details on the columns below, in 'Run parameters'. The path can be either relative or absolute, but cannot contain spaces. 
- `read_files_dir`: A path to a folder containing the read files. Can be relative or absolute. Read files can be either zipped (.fastq.gz) or unzipped (.fastq). Either the containing folder or the filename must begin with the sample identifier. 
  Ex. When processing sample 'A234', the pipeline will find the following: 
  - 'read_files_dir/A234_L001-ds.aaa/other_filename.fastq.gz'
  - 'read_files_dir/A234_L001_R1_001.fastq'
  - 'read_files_dir/other_foldername/A234_L001_R1_001.fastq.gz'


# Run parameters:

Run parameters with information on each sample are provided in a `.csv` file. An example is in `inputs/experiment_info.csv`. Column details:

- `Sample`: Required. The identifier for a sample. It must match the start of the filename or the start of the folder containing the files in the `read_files_dir`, see examples above. This is how sample information is connected to read files. 
- `Description`: Optional. Text describing the sample, it is added to generated charts. 
- `read_type`: Required. How the reads were generated, must be either `restriction` or `fragment`. `restriction` signifies reads were generated via a MmeI restriction enzyme in the transposon end, while `fragment` signifies reads were generated in a random fragmentation process. 
- `transposon_end_sequence`: Required. The sequence of the transposon end (assumed to be right end unless `transposon_end_side` is defined). 
- `transposon_end_side`: Optional. The side of the transposon the end came from: assumed to be the right end unless this is set to `Left`
- `Spacer`: Required. The spacer sequence used in the experiment. 
- `Target fasta file`: Required. A path (either relative or absolute) to a fasta file containing the target molecule's (ex. a genome) sequence. Generates charts of integration sites. 
- `Second target fasta file`: Optional. A path (either relative or absolute) to a fasta file containing a second target's (ex. a second chromosome, a plasmid, etc) sequence. Generates charts of integration sites. 
- `Donor sequence`: Optional. The sequence of the donor molecule. If provided, transposon flanking sequences not mapping to any targets are checked against this sequence, and tallies included in the output logs. 
- `Spike in sequence`: Optional. The sequence of a spike-in. If provided, transposon flanking sequences not mapping to any targets are checked against this sequence, and tallies included in the output logs. 
- `CRISPR Array Sequence`: Optional. The sequence of the CRISPR Array. If provided, transposon flanking sequences not mapping to any targets are checked against this sequence and it's reverse complement, and tallies included in the output logs. 
- `Experiment date`: Optional. Added to chart descriptions. Should be in the format `YYYYMMDD`, if valid this will be used in the output filenames, otherwise the current date is used. 
- `End of protospacer`: Optional. A coordinate number in the target genome used to add a small triangle underneath that point in genome-wide charts. 


# Outputs

The outputs are put in a created `outputs` folder. The `output_log.csv` contains information about each sample run, including read counts and percentages at each filtering step, for matches to donor, spike-ins, etc. The percentages are from reads passing the initial quality filter. The `Experiment date` is copied over from the info_file, if valid, while `Run date` is the date the code was run on. 

In the `outputs/samples` directory, each sample has a CSV with the locations and orientations of all mapped reads. If a second target file file was provided, another CSV is generated with mappings to the second target. These locations are 0-indexed, so will be 1 less than locations as seen in Benchling or other viewers that start at 1. 

A series of plots are generated for each sample in the `outputs/plots` directory, variations of genome-wide mappings with reads bucketed together, and zoomed in plots around the target region. 
