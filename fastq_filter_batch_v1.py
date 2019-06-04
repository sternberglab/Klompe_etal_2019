py_name = 'fastq_filter_batch_v1.py'

###############
# Leo Vo Feb 2019

# Script Function: 3-in-1 script
# 1/ Extracts .gz files from individual subfolders to return fastq files (usually each fastq correspond to 1 lane of 4)
# 2/ Concatenate all separate fastq files for lanes (usually 4) of one sample into 1 main fastq per sample
# 3/ Change all bases with 'low' quality scores to N, and reads with more than half bases with low Q score are filtered.
# 'Low' Qscore is determined by a threshold variable, and function returns a new fastq files with edited filtered reads.
# Keeps track of filtering stats with an excel spreadsheet output

# Script Input(s):
# 1/ One batch info file (CSV) containing numbers 1 to number of files in column 1,
# and Sample IDs (or 'codes') e.g. 'A4715' in column 2. Do not remove header row
# 2/ If running functions separately, note that concat and filtergen fucntions take .fastq files, which are usually
# outputs from previous functions in the script. Functions will output warnings accordingly if input files are not found

# Script Output(s):
# See individual functions. If run all together, returns one filtered fastq file for each SampleID/code/codename,
# and one main excel output for stats keeping purposes.
###############

import os
import csv
import fnmatch
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import xlsxwriter
from pathlib import Path

# change directory based on where input files are
os.chdir('C:\\Users\\Leo Vo\\Desktop\\SEK_NGS\\FASTQ_2019-02-03')

# define run options for the 3 main functionalities of this script
unzip_on = True
concat_on = True
Qfilter_on = True
Qscore_threshold = 20  # bases <= this value are switched to N
excel_out_name = 'fastq_log.xlsx'


# for unzipping .gz files (from individual folders) and moving them into main (current working) folder
def unzip(code):
    for i in os.listdir('.'):
        if fnmatch.fnmatch(i, "{}*".format(code)):
            for x in os.listdir('./{}/'.format(i)):
                name = os.path.splitext(x)[0]
                with gzip.open('./{}/{}'.format(i,x), 'rb') as f_in:
                    with open('./{}'.format(name), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            shutil.rmtree('./{}/'.format(i))  # removes individual folders


# for concatenating FastQ files for indv lanes into one, then moving the indv files into a different folder
def concat(code):
    # search for fastq files containing 5-char 'code at the start', append to collection list
    filenames = []
    for i in os.listdir('.'):
        if fnmatch.fnmatch(i, "{}*.fastq".format(code)):
            filenames.append(i)
    if len(filenames) == 0:
        print('WARNING - NO FILE FOUND FOR {}'.format(code))
        return
    if 0 < len(filenames) < 4:
        print("WARNING FEWER THAN 4 FILES FOUND FOR {}".format(code))
    # concatenate all the files from list end-to-end to a new fastq
    with open('{}_CAT.fastq'.format(code), 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line) # writes new file line by line
            os.rename('./{}'.format(fname), './IndvFastQs/{}'.format(fname))  # move old fastq into subfolder


# initiate excel output containing filtering stats; consider switching to csv in the future
if Qfilter_on:
    excel_out_path = Path(excel_out_name)
    if not excel_out_path.is_file():
        log = xlsxwriter.Workbook(excel_out_name)
        bold = log.add_format({'bold': True})
        logsheet = log.add_worksheet("Qscore Filter Log")
        logsheet.write(0,0, "File Codename", bold)
        logsheet.write(0, 1, "Total Reads", bold)
        logsheet.write(0, 2, "Reads Passing QScore Filter", bold)
        logsheet.write(0, 3, "Percent Reads Filtered", bold)


def tally(file, row, init):  # short function for tallying totals and write to excel output
    total = 0
    for record in SeqIO.parse(file, "fastq"):
        total += 1
    if init == 'init':  # writes to cell for 'pre-filter' total
        logsheet.write(row, 1, total)
    if init == 'fin':  # writes to cell for 'post-filter' total
        logsheet.write(row, 2, total)


def filtergen(file):  # generator function that returns edited reads that pass filter, to write new fastq file
    for record in SeqIO.parse(file, "fastq"):
        # Convert base qualities to Boolean based on Qscore threshold value. Only use reads with >=50% non-N:
        recordqual = [x > Qscore_threshold for x in record.letter_annotations['phred_quality']]  # list of True, False etc
        if float(sum(recordqual)) / float(len(recordqual)) >= .5:  # note that True = 1, False = 0 for summing
            # generates new read sequence where all bases < threshold is switched to 'N'
            seq = "".join([y if x else 'N' for (x, y) in zip(recordqual, record.seq)])
            # create new SeqRecord with edited read sequence
            newrec = SeqRecord(Seq(seq, SingleLetterAlphabet()), id=record.id, name=record.name,
                               description=record.description, letter_annotations=record.letter_annotations)
            yield newrec


def Qfilter(codename, row):  # main command function for filtering
    file_found = 0
    # locates fastq file from codename
    for i in os.listdir('.'):
        if fnmatch.fnmatch(i, "{}*.fastq".format(codename)):
            file = i
            file_found +=1
    if file_found == 1:
        logsheet.write(row, 0, codename)
        tally(file, row, 'init')  # count pre-filter total reads
        # generate input records to write new filtered fastq output file
        SeqIO.write(filtergen(file), '{}_FINAL.fastq'.format(codename), "fastq")
        tally('{}_FINAL.fastq'.format(codename), row, 'fin')  # count post-filter reads in this fastq output
        os.remove('./{}'.format(file))  # remove original fastq CAT input
        return
    elif file_found == 0:
        logsheet.write(row, 0, codename)  # other columns will be left blank
        print("WARNING - No filename found for ", codename)
        return
    else:
        print('WARNING - Multiple Files Found for {}'.format(codename))
        return


# read in the info csv file and begin processing files as a batch
with open('filter_main.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader, None)  # skips header row
    for line in reader:
        row = line[0]
        code = line[1]
        if unzip_on:
            unzip(code)
        if concat_on:
            concat(code)
        if Qfilter_on:
            Qfilter(code, row)

if Qfilter_on:
    log.close()  # close excel output
