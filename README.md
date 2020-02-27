This repo contains the code for interpreting and analyzing Illumina sequencing data. 

The intermediates should be safe to remove (and can be automatically deleted as the script runs, see the parameters), but keep the files in the "outputs" directory so that you can easily graph the results and tweak them to see results. 

Windows requirements:
perl, bowtie2 (2.3.4) and gawk in your path


# in mapping, look first for reads (unique and duplicated), then look for donor or spike matches, then lastly for a CRISPR spacer array sequence from the input csv
