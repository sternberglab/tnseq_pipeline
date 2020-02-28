This repo contains the code for interpreting and analyzing Illumina sequencing data. 

The intermediates should be safe to remove (and can be automatically deleted as the script runs, see the parameters), but keep the files in the "outputs" directory so that you can easily graph the results and tweak them to see results. 

# Getting started:
To run these scripts, you need to be able to run python, perl, bowtie2, and awk. 

### Mac/Linux
For mac/linux, this just means [downloading bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/), any recent version will work, and adding the extracted folder to your PATH environment variable. The others should already be present in your system. 

### Windows
For windows, you'll need to [download Strawberry Perl](http://strawberryperl.com/), download the bowtie2 zip labelled `mingw` [here](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4/) and the `Binaries` zip [here](http://gnuwin32.sourceforge.net/packages/gawk.htm). 

Put the bowtie2 and gawk downloaded folders in a place where you won't need to move them or delete them in the future. Also, for windows machines you should make sure they are located in a full directory path WITHOUT spaces. So if your username contains spaces (C:\Chris Acree\important_downloads\), don't leave the bowtie2 and gawk folders in a subfolder of your user directory. 

Then use the windows search bar to open "Edit the system environment variables", click "Environment Variables". In the following screen, select the "Path" Variable and click "Edit...". Add two new entries to the list, one path to the extracted bowtie2 folder, and one to the bin subfolder of the extracted `gawk` folder (so might be `C:\important_path\gawk-3.1.6-1-bin\bin`). 

## Checking
- To make sure they correctly installed, open Terminal (for Mac) or Command Prompt (for windows) and run each of the following commands:
- python --version
- perl --version
- bowtie2 --version
- awk --version

If any did not print some output showing a version number, it isn't installed. 

To use this code, your python environment also needs to have installed the packages listed in the "Pipfile". Running `pip install -r requirements.txt` from the command line in the project directory will do this. 


# Running this code:
The commonly modified parameters are set in either a CSV for per-sample variables, and in the `parameters.py` file for more general parameters. 

At a minimum, you will need to change the values in parameters.py:
- Set `working_dir` to a directory that will be modified and added to by the script. It should have one subfolder `raw` containing either the zipped or raw fastq files from your Illumina runs. 
- Set `genome_path` to point to your genome fasta file to run alignment against. 
- Set `info_file` to point to the information csv file with meta information about each sample run. There is an example of the required format for this csv in `inputs/experiment_info`, you can choose to modify that file and use it directly or make a new csv with the same columns and have the `info_file` path point to it. 

The information CSV file must exist in the specific format with information about each sample being analyzed. Missing values are likely to cause failures. 


# Outputs
The outputs are put in a created `outputs` folder under the working directory set in parameters.py The `output_log` contains information about each run from raw files, and for each run, a subfolder of sample and Q score is created containing
- A csv with counts of unique reads within the genome
- A fasta containing those unique reads
- An excel file with details about the transposition distance mapping reads
- 3 binned histogram plots of alignment results across the genome, one with raw counts, one with normalized percentages, and one with normalized percentages zoomed in to show lower-level read locations
- 2 localized plots around the target location, showing specific read counts around the target, one as an overlapped histogram and the other as separate plots for RL vs LR reads. 