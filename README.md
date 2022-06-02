This repository contains the code for interpreting and analyzing Illumina sequencing data for mapping transposon insertion events. 

An "outputs" folder is created in this directory with generated plots and CSV files of reads for each processed sample. 


# Getting started:
To run these scripts, you need to be able to run python, perl, bowtie2, and awk. 

### Mac/Linux
For mac/linux, this just means [downloading bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/), any recent version will work, and adding the extracted folder to your PATH environment variable. The others should already be present in your system. 

### Windows
For windows, you'll need to [download Strawberry Perl](http://strawberryperl.com/), download the bowtie2 zip labelled `mingw` [here](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4/) and the `Binaries` zip [here](http://gnuwin32.sourceforge.net/packages/gawk.htm). 

Put the bowtie2 and gawk downloaded folders in a place where you won't need to move them or delete them in the future. Also, for windows machines you should make sure they are located in a full directory path WITHOUT spaces. So if your username contains spaces (eg. C:\Chris Acree\important_downloads\), don't leave the bowtie2 and gawk folders in a subfolder of your user directory. 

Then use the windows search bar to open "Edit the system environment variables", click "Environment Variables". In the following screen, select the "Path" Variable and click "Edit...". Add two new entries to the list, one path to the extracted bowtie2 folder, and one to the bin subfolder of the extracted `gawk` folder (so might be `C:\important_path\gawk-3.1.6-1-bin\bin`). 

## Checking
- To make sure they correctly installed, open Terminal (for Mac) or Command Prompt (for windows) and run each of the following commands:
- python --version
- perl --version
- bowtie2 --version
- awk --version

If any did not print some output showing a version number, it isn't installed properly. 

To use this code, your python environment also needs to have installed the packages listed in the "Pipfile". Running `pip install -r requirements.txt` from the command line in the project directory will do this. 


# Running this code:
The commonly modified parameters are set in a CSV for per-sample variables, and in the `parameters.py` file for more general but less frequently changed parameters. 

At a minimum, you will need to change the values in `parameters.py`:
- Set `working_dir` to a directory containing a `raw` subfolder with either the zipped or raw fastq files from your Illumina runs. 
- Set `info_file` to point to the information csv file with meta information about each sample run. There is an example of the required format for this csv in `inputs/experiment_info`, you can choose to modify that file and use it directly or make a new csv with the same columns and have the `info_file` path point to it. 

The information CSV file must exist in the specific format with information about each sample being analyzed. For each sample, the minimum required fields are `Sample`, `read_type` (either "restriction" or "fragment"), `tn_end_sequence`, `Spacer`, and `Genome` (a full path to the fasta file). 

Adding `CRISPR Array Sequence`, `Donor sequence` and `Spike in sequence` will generate additional statistics and, if plasmid fasta file is also provided, plots. 

`End of protospacer` attaches small triangles to the genome wide plots at the given location. 

If you're parameters and input file are ready, run the `main.py` file. 


# Outputs
The outputs are put in a created `outputs` folder. The `output_log` contains information about each sample run from raw files like the reads passing each step. In the `outputs/samples` directory, each sample has a CSV with the locations and orientations of all mapped reads. If plasmid file was provided, another CSV is generated with mappings to the plasmid. 

A series of plots are generated for each sample in the `outputs/plots` directory, variations of genome-wide mappings with reads bucketed together, and zoomed in plots around the target region. 
