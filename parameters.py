# 1. FOR LOADING NEW DATA

# The directory where your raw files are held (can be zip files)
# This can be either relative or absolute path
rawFilesDir = "../raw_files"

# The Qscore_threshold to use
# This must be defined (use -1 if do not want to filter)
Qscore_threshold = 20


# 2. FOR FINGERPRINTING

# The R flank sequence to look for
flank_sequence = "TGTTGGAACAACCAT"

# The length of the fingerprint to match (set by restriction enzyme)
fingerprint_length = 17


# 3. 