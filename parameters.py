# General settings

# This will delete intermediate files
# It will make full analysis from scratch take longer
# but free up disk space
delete_intermediates = True

# Path to an info csv file with various information, see the example file for details
info_file = "../rename/input_selective.csv"

# 1. FOR LOADING NEW DATA

# The directory where your raw Illumina files are held, 
# as well as where intermediaries and outputs will be saved
# Illumina files should be placed in a /raw subfolder
# and can be zip files or unzipped
# This can be either relative or absolute path
working_dir = "../rename"

# The Qscore_threshold to use
# This must be defined (use -1 if do not want to filter)
Qscore_threshold = 20


# 2. FOR FINGERPRINTING

# The R flank sequence to look for
flank_sequence = "TGTTGGAACAACCAT"

# The length of the fingerprint to match (set by restriction enzyme)
fingerprint_length = 17


# 3. FOR MAPPING

# The number of bases duplicated by transposon insertion
transposon_site_duplication_length = 5

# 4. GRAPHS
# For all graphs, use the following filetype and dpi
plots_filetype = 'svg'
plots_dpi = 600
fig_size_inches = [8.9, 3.1]

# Bin size for genome-wide plots in base pairs
genome_bin_size = 5000

# Percent of reads at which to cap the y-axis on the zoomed in histogram
# to show low-level reads
low_reads_cap_percent = 0.50

# Run parameters for transposition distance histogram
query_length = 500
on_target_window = 100