# General settings

#### REQUIRED ####

'''
Path to the csv file with information about the samples to run.
See an example in './inputs/experiment_info.csv', and descriptions
of columns are in the README. 
Paths can be either relative or absolute, but cannot contain spaces. 
'''
info_file = "./experiment_info.csv"

'''
Path to a folder containing the read files. Can be relative or absolute. 
Read files can be either lzipped (.fastq.gz) or unzipped (.fastq)
Either the containing folder or the filename must begin with
the sample identifier. 
Ex. When processing sample 'A234', the pipeline will find the following: 
- 'read_files_dir/A234_L001-ds.aaa/other_filename.fastq.gz'
- 'read_files_dir/A234_L001_R1_001.fastq'
- 'read_files_dir/other_foldername/A234_L001_R1_001.fastq.gz'
'''
read_files_dir = '../ngs_run_220711/'


#### DEBUGGING SETTINGS ####

'''
Intermediate files are created at different steps of processing. 
Setting this to false will keep these files for debugging or looking
deeper at some filtered reads. 
Normally, these files are deleted to save space. 
'''
delete_intermediates = True


#### SENSITIVITY AND SYSTEM SETTINGS ####

'''
The Qscore_threshold to use. This must be defined. 
Use -1 to not filter any reads. 
Reads with <50% bases passing threshold are filtered. 
'''
Qscore_threshold = 20

# The length of the sequence to match to the target files
transposon_end_flanking_sequence_length = 17

# The number of bases duplicated by transposon insertion
target_site_duplication_length = 5


#### GRAPHING SETTINGS ####

# General matplotlib settings
plots_filetype = 'svg'
plots_dpi = 600
fig_size_inches = [8.9, 3.1]

# Bin size for genome-wide plots, in base pairs
genome_bin_size = 5000

# Percent of reads at which to cap the y-axis on the zoomed in histogram
# to show low-level reads
low_reads_cap_percent = 0.50

# On-target calculated window in the transposition distance histogram
on_target_window = 200

# Distance from the end of the protospacer to the target
dist_to_target = 50
