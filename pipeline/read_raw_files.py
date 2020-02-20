from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from .utils import inter_path

def filtergen(file, threshold):  # generator function that returns edited reads that pass filter, to write new fastq file
    for record in SeqIO.parse(file, "fastq"):
        # Convert base qualities to Boolean based on Qscore threshold value. Only use reads with >=50% non-N:
        recordqual = [x > threshold for x in record.letter_annotations['phred_quality']]  # list of True, False etc
        if float(sum(recordqual)) / float(len(recordqual)) >= .5:  # note that True = 1, False = 0 for summing
            # create new SeqRecord with edited read sequence
            newrec = SeqRecord(record.seq, id=record.id, name=record.name,
                               description=record.description, letter_annotations=record.letter_annotations)
            yield newrec

def process_files(input_filenames, filtered_path):
    total_records = 0
    concat_path = inter_path('{}_CONCAT.fastq'.format(code))
    with open(inter_path(concat_path), 'w') as outfile:
        for file in input_filenames:
            with open(file, 'r') as fastq_file:
                for line in fastq_file:
                    outfile.write(line)
                    if line.startswith('@'):
                        total_records += 1

        with open(filtered_path, 'w+') as outfile:
            filtered_records = SeqIO.write(filtergen(concat_path, Qscore_threshold), outfile, "fastq")

    os.remove(concat_path)
    return (total_records, filtered_records)