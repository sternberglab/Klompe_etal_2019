py_name = 'trans_dist_dual_anchor_v1.py'

# Leo Vo Feb 2019
# For measuring distance between end of protospacer to TN7 flanks for NGS reads from a dual-anchored primers method.

import os
import fnmatch
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
import xlsxwriter
from pathlib import Path
import datetime as dt

# change if necessary
os.chdir('C:\\Users\\Leo Vo\\Desktop\\SEK_NGS\\Distance')

# change if necessary
exp_date = "FEB0219"
user = "Leo Vo"

tn7_rflank = 'TGTTGATACAACCATAAAATGATAATTACACCCATAAATTGATAATTATC'
tn7_lflank = 'TGTTGATGCAACCATAAAGTGATATTTAATAATTATTTATAATCAGCAAC'
donor_R = str('gataacaatttcacacaggaaacagctatgaccatgattacgccaagctt').upper()
donor_R = Seq(donor_R, SingleLetterAlphabet())
donor_L = str('cccagtcacgacgttgtaaaacgacggccagtgaattcgagctcggtacc').upper()
donor_L = Seq(donor_L, SingleLetterAlphabet())


#initializes excel output file
excel_out_path = Path('{}_transposition_distance_dual_anchor.xlsx'.format(exp_date))
if not excel_out_path.is_file():
    log = xlsxwriter.Workbook('{}transposition_distance_dual_anchor.xlsx'.format(exp_date))
    bold = log.add_format({'bold': True})
    upsizebold = log.add_format()
    upsizebold.set_font_size(18)
    upsizebold.set_bold()
    percentage_format = log.add_format()
    percentage_format.set_num_format('0.00%')
    deci3_format = log.add_format()
    deci3_format.set_num_format('0.000')
#make landing info sheet containing all information shared by sheets in the excel file
infosheet = log.add_worksheet("General Info")
infosheet.write(0, 0, "NGS Dual Anchored Analysis of Transposition Distance ", upsizebold)
infosheet.set_column(0, 0, 35)
infosheet.write(2, 0, "NextSeq/NGS Data Collection Date", bold)
infosheet.write(2, 1, exp_date)

infosheet.write(3, 0, "Python Data Analysis Date", bold)
infosheet.write(3, 1, str(dt.datetime.now())[0:10])

infosheet.write(4, 0, "Username/Initials", bold)
infosheet.write(4, 1, user)

infosheet.write(5, 0, "Python Code Used", bold)
infosheet.write(5, 1, py_name)

tn7_rflank = 'TGTTGATACAACCATAAAATGATAATTACACCCATAAATTGATAATTATC'
tn7_lflank = 'TGTTGATGCAACCATAAAGTGATATTTAATAATTATTTATAATCAGCAAC'

def distance_find(record, read_flank, flank_length, geno_length, opposite, refseq_primer, geno_start):
    # locating transposon flank
    flank_start = -1
    if read_flank == 'R':
        if record.seq.find(tn7_rflank[0:flank_length]) >= 0:
            flank_start = record.seq.find(tn7_rflank[0:flank_length])
        if record.seq.find(tn7_lflank[0:flank_length]) >= 0:
            flank_start = record.seq.find(tn7_lflank[0:flank_length])
            print("WARNING - Found L flank instead of R")
    if read_flank == 'L':
        if record.seq.find(tn7_lflank[0:flank_length]) >= 0:
            flank_start = record.seq.find(tn7_lflank[0:flank_length])
        if record.seq.find(tn7_rflank[0:flank_length]) >= 0:
            flank_start = record.seq.find(tn7_rflank[0:flank_length])
            print("WARNING - Found R flank instead of L")

    # map adjacent genomic sequence to ref genome
    # important to RESET these Boolean switches right here otherwise they read forward into other reads
    donor = False
    trans_dist = -1
    N_found = False
    flank_found = False

    if flank_start > 0:
        flank_found = True
        genome_map_seq = record.seq[flank_start - geno_length:flank_start].upper()
        if refseq_primer.find(genome_map_seq) >= 0:
            # returns distance from end of protospacer
            if not opposite:
                trans_dist = refseq_primer.find(genome_map_seq) + geno_length - geno_start
            if opposite:
                trans_dist = len(refseq_primer) - refseq_primer.find(genome_map_seq) - geno_length + 5 - geno_start
        if donor_L.find(genome_map_seq) >= 0 or donor_R.find(genome_map_seq) >= 0:
            donor = True
        if genome_map_seq.find('N') >= 0:
            N_found = True
    else:
        flank_found = False
    return trans_dist, flank_found, record.seq, donor, N_found


def analyze_file():
    # inputs from user
    print("NOTE - All user inputs case insensitive")
    filename = input("Enter FastQ Sample ID (e.g. A804): ").upper()
    for i in os.listdir('.'):
        if fnmatch.fnmatch(i, "{}*.fastq".format(filename)):
            file = i
    description = input("Enter sample description: ").upper()
    refseq = input("Enter 1kb E.coli reference genomic sequence: ").upper()
    refseq = Seq(refseq, SingleLetterAlphabet())
    pam = input("Enter the PAM sequence: ").upper()
    spacer = input("Enter the protospacer sequence ").upper()
    read_flank = input("Enter which Tn7 flank is included in the i5 read ('R' for right, 'L' for left): ").upper()
    t_location = input("Enter whether protospacer is on the same forward strand of the reference genomic"
                       "sequence ('FW') or on the reverse strand ('RV'): ").upper()
    p_location = input("Enter whether genomic i5 primer is on the same forward strand of the reference genomic"
                       "sequence ('FW') or on the reverse strand ('RV'): ").upper()
    geno_p = input("Enter oSL number for genomic primer WITHOUT oSL part: ")
    tran_p = input("Enter oSL number for transposon primer WITHOUT oSL part: ")
    flank_length = int(input("Enter length of Tn7 flank used for mapping reads (up to 40): "))
    geno_length = int(
        input("Enter length of genomic sequence (adjacent to Tn7 flank) used for mapping reads (up to 40): "))

    error_log = []  # stores errors for excel output later

    # checks PAM and protospacer and set up coordinate system starting at the end 3' of protospacer
    if t_location == 'FW':
        refseq_spacer = refseq
    if t_location == 'RV':
        refseq_spacer = refseq.reverse_complement()
    if p_location == 'FW':
        refseq_primer = refseq
    if p_location == 'RV':
        refseq_primer = refseq.reverse_complement()

    if t_location == p_location:
        opposite = False
    else:
        opposite = True

    if refseq_spacer.find(spacer) >= 0:
        geno_start = refseq_spacer.find(spacer) + len(spacer)
        pam_start = refseq_spacer.find(spacer) - 2
        # confirm PAM
        if refseq_spacer[pam_start:(pam_start + 2)] != pam:
            print("WARNING - Incorrect PAM")
            error_log.append("WARNING - Incorrect PAM")
    else:
        print("ERROR - Spacer not found!")
        return

    out_tally = [0]*len(refseq)
    seq_list = ['X']*len(refseq)
    non_maps = 0
    total = 0
    flank_count = 0
    donor_count = 0
    n_count = 0

    for record in SeqIO.parse(file, "fastq"):
        total += 1
        trans_dist, flank_found, seq, donor, N_found = distance_find(record, read_flank, flank_length, geno_length,
                                                            opposite, refseq_primer, geno_start)
        if trans_dist >= 0:
            if out_tally[trans_dist] == 0:
                out_tally[trans_dist] += 1
                seq_list[trans_dist] = str(seq)
            else:
                out_tally[trans_dist] += 1
        else:
            non_maps += 1

        if flank_found:
            flank_count += 1
        if donor:
            donor_count += 1
        if N_found:
            n_count += 1
    #sets up excel output sheet
    logsheet = log.add_worksheet(filename)

    logsheet.set_column(0, 0, 32)
    logsheet.write(0, 0, "Sample ID", bold)
    logsheet.write(0, 1, filename)

    logsheet.write(1, 0, "Description", bold)
    logsheet.write(1, 1, description)

    logsheet.write(2, 0, "Target Location", bold)
    if t_location == 'FW':
        logsheet.write(2, 1, "5' of Integration Site")
    else:
        logsheet.write(2, 1, "3' of Integration Site (RevCom)")

    logsheet.write(3, 0, "Primer Location", bold)
    if p_location == 'FW':
        logsheet.write(3, 1, "5' of Integration Site")
    else:
        logsheet.write(3, 1, "3' of Integration Site (RevCom)")

    logsheet.write(4, 0, "Transposon End Being Queried", bold)
    logsheet.write(4, 1, read_flank)

    logsheet.write(5, 0, "Genomic Primer ID", bold)
    logsheet.write(5, 1, "oSL{}".format(geno_p))

    logsheet.write(6, 0, "Transposon Primer ID", bold)
    logsheet.write(6, 1, "oSL{}".format(tran_p))

    logsheet.write(7, 0, "Full Genomic Query Sequence", bold)
    logsheet.write(7, 1, str(refseq))

    logsheet.write(8, 0, "PAM Sequence", bold)
    logsheet.write(8, 1, pam)

    logsheet.write(9, 0, "Protospacer Sequence", bold)
    logsheet.write(9, 1, spacer)

    logsheet.write(10, 0, "Transposon End Filter Length", bold)
    logsheet.write(10, 1, flank_length)

    logsheet.write(11, 0, "Genomic Sequence Query Length", bold)
    logsheet.write(11, 1, geno_length)

    logsheet.write(12, 0, "Total Reads", bold)
    logsheet.write(12, 1, total)

    logsheet.write(13, 0, "Fully Mapped Reads", bold)
    logsheet.write(13, 1, sum(out_tally))
    logsheet.write(13, 2, sum(out_tally)/total, percentage_format)

    logsheet.write(14, 0, "Unmapped Reads", bold)
    logsheet.write(14, 1, non_maps)
    logsheet.write(14, 2, non_maps/total, percentage_format)

    logsheet.write(15, 0, "Unmapped Reads with Mapped Transposon End", bold)
    logsheet.write(15, 1, flank_count)
    logsheet.write(15, 2, flank_count/total, percentage_format)

    logsheet.write(16, 0, "Genomic Queries Containing N(s)", bold)
    logsheet.write(16, 1, n_count)
    logsheet.write(16, 2, n_count / total, percentage_format)

    logsheet.write(17, 0, "Reads Mapping to Donor Plasmid", bold)
    logsheet.write(17, 1, donor_count)
    logsheet.write(17, 2, donor_count / total, percentage_format)

    logsheet.write(18, 1, "Protospacer-Transposon Distance", bold)
    logsheet.write(18, 0, "Genomic Base", bold)
    logsheet.write(18, 2, "Number of Reads", bold)
    for i in range(0, len(refseq)):
        logsheet.write(i + 19, 1, i)
    for i in range(0, len(refseq)):
        logsheet.write(i + 19, 0, refseq[i])
    for i in range(0, len(out_tally)):
        logsheet.write(i + 19, 2, out_tally[i])

    logsheet.write(18, 3, "% of Total Mapped Reads", bold)
    for i in range(0, len(out_tally)):
        logsheet.write(i + 19, 3, out_tally[i]/sum(out_tally), percentage_format)

    logsheet.write(18, 4, "Normalized Read Count", bold)
    for i in range(0, len(out_tally)):
        logsheet.write(i + 19, 4, out_tally[i] / max(out_tally), deci3_format)

    logsheet.write(18, 5, "Example Read Sequence", bold)
    for i in range(0, len(seq_list)):
        logsheet.write(i + 19, 5, seq_list[i])
    # writes error_log to a txt file
    #error_file = open("{}_error_log.txt".format(filename), "w")
    #for i in error_log:
    #    error_file.write(i)
    #error_file.close()


analyze_file()

log.close()