py_name = 'site_trans_mmei_v1B_batch.py'

###############
# Leo Vo Feb 2019

# Script Function:
# Process filtered reads outputted from Geneious and map to reference genome. Provide locations of reads based on
# 3' edge of each mapping sequence (each read is only tallied once at the genomic position corresponding to 3' edge)
# V1B: Supports duplicated reads

# Script Input(s):
# 1/ One batch info file (CSV) containing (at least) a list of Sample IDs/Codes (e.g A4750)
# 2/ Corresponding fastq files outputted from Geneious, which contain reads of slightly varying length which have
# been confirmed to map to the reference genome. These files have corresponding Sample ID in their names.
# 3/ A genome.fasta file containing the ecoli refseq
# 4/ Make sure subfolder 'Trans_Sites\' exists

# Script Output(s):
# One .csv file for each input fastq, tallying number of reads corresponding to each genomic coordinate.
###############

import os
import fnmatch
from Bio import SeqIO
import csv

# change directory based on where input files are
os.chdir('C:\\Users\\Leo Vo\\Desktop\\SEK_NGS\\190217_Integration_site_analysis_SEK')

# some run parameters
exp_date = "190202"
map_length = 17  # length of final trimmed sequence used for mapping. i.e. 17 for Tn7, 14 for Himar
TSD = 5  # TSD length for TSD correction (likely not necessary due to subsequent binning of data)

# read in BL21DE3 RefSeq
for record in SeqIO.parse("genome.fasta", "fasta"):
    genome = record.seq.upper()  # remember to convert to upper case
    genome_rc = genome.reverse_complement()
len_genome = len(genome)
# add bases to the end from front to correct for sequences at the linear junctions
genome = genome + genome[:map_length-1]
genome_rc = genome_rc + genome_rc[:map_length-1]

# function to find all occurences of each read
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)  # use start += 1 to find overlapping matches


def site_find(code, filename, genome, genome_rc, len_genome, map_length):  # main analysis function
    total = 0  # counts total reads in fastq file
    fw = 0  # counts total forward reads, not currently used for anything
    rv = 0  # counts all reverse reads, not currently used for anything
    tally = [0]*len_genome  # list for tallying locations of transposition
    for record in SeqIO.parse(filename, "fastq"):  # loop through all reads in fastq file
        total += 1
        rev_seq = record.seq.reverse_complement()  # reverse the read (artifact of Geneious output)
        new_seq = rev_seq[len(rev_seq) - map_length:].upper()  # trim to last map_length base pair

        # map to genome and tally
        locations = list(find_all(genome, new_seq))  # look in forward strand first
        if len(locations) > 0:
            fw += len(locations)
            for i in locations:
                site = i + map_length
                if site > len_genome:  # correct for junction indexing
                    site = site - len_genome
                tally[site - 1] += 1  # -1 to correct for indexing

        locations_rc = list(find_all(genome_rc, new_seq))  # then look in reverse strand
        if len(locations_rc) > 0:
            rv += len(locations_rc)
            for m in locations_rc:
                site_rc = len_genome - m - map_length + TSD  # convert to fwd strand coordinates
                if site_rc < 0:
                    site_rc = site_rc + len_genome
                tally[site_rc - 1] += 1  # -1 to correct for indexing

    # output to .csv file (2 columns - Position & Coverage)
    # make sure subfolder 'Trans_Sites' exists
    with open('Trans_Sites\\{}_Trans_Sites.csv'.format(code), 'w', newline='') as csv_out:
        writer = csv.writer(csv_out)
        writer.writerow(['Total Reads', total])  # include row with total reads
        writer.writerow(['Position', 'Coverage'])  # include header row
        for i in range(0, len(tally)):
            writer.writerow([i+1, tally[i]])
    csv_out.close()


# use batch .csv input file to run all files back-to-back
with open('input_spacers_test.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader, None)
    for row in reader:
        code = row[0]
        filename = 'none'
        for i in os.listdir('.'):
            if fnmatch.fnmatch(i, "*{}*.fastq".format(code)):  # search for file with matching Sample ID in cwd
                filename = i
                break
        if filename != 'none':
            # run main site_find function with the relevant input variables
            site_find(code, filename, genome, genome_rc, len_genome, map_length)
        else:
            print("WARNING - File Not Found For {}".format(code))
