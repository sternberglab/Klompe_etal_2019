py_name = 'trans_dist_mmei_v1_batch.py'

# Leo Vo, Feb 2019
# For measuring distance between end of protospacer to TN7 flanks for NGS reads from an MmeI NGS library
# Reads in FastQ files outputed from GENEIOUS Prime Software, which have transposon flank sequences removed
# and are reverse complemented from the original Illumina fastq output.
# Uses a master csv to run multiple files in a batch
# Outputs: The script can output multiple things depending on the output options in the code below:
#          - An Excel file containing all relevant information
#          - Plots of transposition distances; either RL and LR plotted separately, or overlapped together

import os
import fnmatch
from Bio import SeqIO
import xlsxwriter
from pathlib import Path
import datetime as dt
import csv
import heapq
import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'  # important so that text stays as characters in the svg output

# change directory based on where input files are
os.chdir('C:\\Users\\Leo Vo\\Desktop\\SEK_NGS\\190217_Integration_site_analysis_SEK')

# change if necessary
exp_date = "190202"
user = "Leo Vo"

# output options
excel = False  # outputs an excel spreadsheet summary
excel_name = '{}_transposition_distance_mmei_duplicates.xlsx'.format(exp_date) ## Change output excel file name here
plot = True  # outputs plots of integration around primary site
plot_overlap = True  # overlaps the distributions of the opposite orientations **Will not plot if plot = False above

#reads in BL21DE3 RefSeq
for record in SeqIO.parse("genome.fasta", "fasta"):
    genome = record.seq.upper()  # remember to convert to upper case
    genome_rc = genome.reverse_complement()

if excel:
    # initializes excel output file, define formats used below
    excel_out_path = Path(excel_name)
    if not excel_out_path.is_file():
        log = xlsxwriter.Workbook(excel_name)
        bold = log.add_format({'bold': True})
        upsizebold = log.add_format()
        upsizebold.set_font_size(16)
        upsizebold.set_bold()
        percentage_format = log.add_format()
        percentage_format.set_num_format('0.00%')
        deci3_format = log.add_format()
        deci3_format.set_num_format('0.000')
        red_font = log.add_format({'font_color': 'red'})
        blue_font = log.add_format({'font_color': 'blue'})
        red_percent = log.add_format()
        red_percent.set_num_format('0.00%')
        red_percent.set_font_color('red')
        blue_percent = log.add_format()
        blue_percent.set_num_format('0.00%')
        blue_percent.set_font_color('blue')
        red_deci3 = log.add_format()
        red_deci3.set_num_format('0.000')
        red_deci3.set_font_color('red')
        blue_deci3 = log.add_format()
        blue_deci3.set_num_format('0.000')
        blue_deci3.set_font_color('blue')

    # make landing info sheet containing all information shared by sheets in the excel file
    infosheet = log.add_worksheet("General Info")
    infosheet.write(0, 0, "NGS MmeI Analysis of Transposition Distance ", upsizebold)
    infosheet.set_column(0, 0, 35)
    infosheet.write(2, 0, "NextSeq/NGS Data Collection Date", bold)
    infosheet.write(2, 1, exp_date)

    infosheet.write(3, 0, "Python Data Analysis Date", bold)
    infosheet.write(3, 1, str(dt.datetime.now())[0:10])

    infosheet.write(4, 0, "Username/Initials", bold)
    infosheet.write(4, 1, user)

    infosheet.write(5, 0, "Python Code Used", bold)
    infosheet.write(5, 1, py_name)

    infosheet.write(6, 0, "Notes", bold)
    infosheet.write(7, 0, '"Query Window" - larger window for query distances (not the same as on-target window'
                          ', same strand as the BL21 RefSeq')
    infosheet.write(8, 0, '"Genomic Base" - Base directly 5\' of integration site, on the SAME strand as protospacer')
    infosheet.write(9, 0, '"Example Reads" - Reads here are trimmed to 17bp and the reverse complement of reads in'
                          'the Geneious fastq output', )

# read in csv file containing sample codes, info and spacers.
# locate fastq file based on code

# for each file, determine spacer orientation then switch genome to revcom if needed (refseq)
# match spacer with refseq
# That + 32 = spacer_end. From refseq take window 400 before and 400 after spacer_end as refseq_query


# for each fastq take reads, reverse and trim them down to the last 17 bp. Then map to refseq_query or revcom of that
# if it maps to location on query, then it's RL, and trans_dist is location + 17 - spacer_end
# if it maps to location on query revcom, it's LR, then trans_dist is len refseq - location -17 + 5 - spacer_end
# add these trans_dist into separate out_tally lists for RL and LR, and also into a common out_tally
# determine most freq trans_dist from common out_tally - main_site
# set on_target as the window -50 and +50 from trans_site.
# go through the separate out tally lists and remove any trans_dist not within on_target
# sum of number of items from these two lists versus total of reads in the fastq is on-target freq
# ratio of the 2 lists is the orientation bias
# output the separate trans dists into excel output in different colors, and also a common column
# use matplotlib to plot and save separately graphs for the 2 orientations


def dist_find(code, psl, description, filename, direction, refseq, spacer, date, excel, plot, plot_overlap):  # main analysis function
    # map spacer to refseq and determine query window
    query_length = 500
    if refseq.find(spacer) >= 0:
        spacer_end = refseq.find(spacer) + 32
        query = refseq[spacer_end-90:spacer_end+query_length]  # '-99' accounts for the -50bp of the on-target later
        query_rc = query.reverse_complement()
        spacer_end = 90  # resets spacer end index to middle of query window (no longer using full refseq)
    else:
        print("ERROR - Spacer not found within RefSeq")
        return
    total = 0  # counts total reads in fastq file
    out_list_all = []  # list holding common trans_dist
    out_list_rl = []  # list holding indv RL trans_dist values
    out_list_lr = []  # list holding indv LR trans_dist values
    # these lists are longer than query_length to hold negative values of trans_dist
    # the output excel cuts of list using query_length so those values will not show
    out_tally_rl = [0] * (query_length+spacer_end+20)  # list tallying freq of tran_dist for RL
    out_tally_lr = [0] * (query_length+spacer_end+20)  # list tallying freq of tran_dist for LR
    example_reads_rl = ['X'] * (query_length+spacer_end+20)  # to hold example reads mapping to each trans_dist for RL
    example_reads_lr = ['X'] * (query_length+spacer_end+20)  # to hold example reads mapping to each trans_dist for LR
    for record in SeqIO.parse(filename, "fastq"):  # loop through all reads in fastq file
        total += 1
        rev_seq = record.seq.reverse_complement()
        new_seq = rev_seq[len(rev_seq)-17:]  # trim to last 17 base pair
        if query.find(new_seq) >= 0:  # corresponds to RL
            trans_dist = query.find(new_seq) + 17 - spacer_end  # distance in bp from end of protospacer
            out_list_all.append(trans_dist)  # append to holding lists for processing later
            out_list_rl.append(trans_dist)
            if out_tally_rl[trans_dist] == 0:  # add read to example list if this is the first occurrence
                example_reads_rl[trans_dist] = new_seq
            out_tally_rl[trans_dist] += 1  # count into tally list
        elif query_rc.find(new_seq) >= 0:  # corresponds to LR
            trans_dist = len(query) - query_rc.find(new_seq) - 17 + 5 - spacer_end  # dist in bp from end of protospacer
            out_list_all.append(trans_dist)  # append to tally lists for processing later
            out_list_lr.append(trans_dist)
            if out_tally_lr[trans_dist] == 0:  # add read to example list if this is the first occurrence
                example_reads_lr[trans_dist] = new_seq
            out_tally_lr[trans_dist] += 1  # count into tally list

    #  determine most frequent trans_dist
    out_tally_all = [0] * query_length
    for i in range(0, len(out_tally_all)):
        out_tally_all[i] = out_tally_rl[i] + out_tally_lr[i]
    for x, y in enumerate(out_tally_all):
        if y == max(out_tally_all):
            main_site = x + spacer_end  # remember to convert dist to site of integration

    # define on target window
    on_target_lower = main_site - 50
    on_target_upper = main_site + 50
    # move any trans_dist within this window into a final holding list and clears old holding list
    final_list_rl = []
    for dist in out_list_rl:
        if on_target_lower <= (dist + spacer_end) <= on_target_upper:  # convert dist to site of integration
            final_list_rl.append(dist)

    final_list_lr = []
    for dist in out_list_lr:
        if on_target_lower <= (dist + spacer_end) <= on_target_upper:  # convert dist to site of integration
            final_list_lr.append(dist)

    # determine on target frequency
    on_target_total = len(final_list_rl) + len(final_list_lr)
    off_target = total - on_target_total

    # determine top 3 most common trans_dist for highlight box
    # for combined RL and LR
    indices = []  # for zipping with out_tally lists
    for i in range(0, query_length):
        indices.append(i)
    top_3 = heapq.nlargest(3, zip(out_tally_all, indices))  # exists as a list of smaller 2-item lists

    if excel:
        #set up excel output sheet
        logsheet = log.add_worksheet(code)

        logsheet.set_column(0, 0, 24)
        logsheet.set_column(1, 1, 20)
        logsheet.set_column(2, 6, 17)
        logsheet.set_column(7, 12, 19)
        logsheet.write(3, 3, " ")  # for clearer aesthetic

        logsheet.write(0, 0, "Sample ID", bold)
        logsheet.write(0, 1, code)

        logsheet.write(1, 0, "Description", bold)
        logsheet.write(1, 1, description)

        logsheet.write(2, 0, "Target Location", bold)
        if direction == 'fw':
            logsheet.write(2, 1, "5' of Integration Site")
        else:
            logsheet.write(2, 1, "3' of Integration Site (RevCom)")

        logsheet.write(3, 0, "Query Window", bold)
        if direction == 'fw':
            logsheet.write(3, 1, str(query))
        if direction == 'rv':
            logsheet.write(3, 1, str(query_rc))

        logsheet.write(4, 0, "Plasmid encoding gRNA", bold)
        logsheet.write(4, 1, psl)

        logsheet.write(5, 0, "Protospacer", bold)
        logsheet.write(5, 1, spacer)

        logsheet.write(6, 0, "Total Reads", bold)
        logsheet.write(6, 1, total)

        logsheet.write(7, 0, "On Target Reads", bold)
        logsheet.write(7, 1, on_target_total)
        logsheet.write(7, 2, on_target_total / total, percentage_format)

        logsheet.write(8, 0, "Off Target Reads", bold)
        logsheet.write(8, 1, off_target)
        logsheet.write(8, 2, off_target / total, percentage_format)

        logsheet.write(9, 0, "On Target Reads in RL Orientation", bold)
        logsheet.write(9, 1, len(final_list_rl))
        logsheet.write(9, 2, len(final_list_rl) / total, percentage_format)

        logsheet.write(10, 0, "On Target Reads in LR Orientation", bold)
        logsheet.write(10, 1, len(final_list_lr))
        logsheet.write(10, 2, len(final_list_lr) / total, percentage_format)

        logsheet.write(11, 1, "Protospacer-Transposon Distance", bold)
        for i in range(0, query_length):
            logsheet.write(i + 12, 1, i)

        logsheet.write(11, 0, "Genomic Base", bold)
        for i in range(-1, query_length - 1):
            logsheet.write(i + 13, 0, query[i+spacer_end])  # shift back 1 to get the base right before transposition

        logsheet.write(11, 2, "Number of Reads (RL)", bold)
        for i in range(0, query_length):
            logsheet.write(i + 12, 2, out_tally_rl[i], red_font)
        logsheet.write(11, 3, "% of Total Reads (RL)", bold)
        for i in range(0, query_length):
            logsheet.write(i + 12, 3, out_tally_rl[i]/total, red_percent)
        logsheet.write(11, 4, "Normalized Read Count (RL)", bold)
        if max(out_tally_rl) > 0:
            for i in range(0, query_length):
                logsheet.write(i + 12, 4, out_tally_rl[i] / max(out_tally_rl), red_deci3)

        logsheet.write(11, 5, "Number of Reads (LR)", bold)
        for i in range(0, query_length):
            logsheet.write(i + 12, 5, out_tally_lr[i], blue_font)
        logsheet.write(11, 6, "% of Total Reads (LR)", bold)
        for i in range(0, query_length):
            logsheet.write(i + 12, 6, out_tally_lr[i] / total, blue_percent)
        logsheet.write(11, 7, "Normalized Read Count (LR)", bold)
        if max(out_tally_lr) > 0:
            for i in range(0, query_length):
                logsheet.write(i + 12, 7, out_tally_lr[i] / max(out_tally_lr), blue_deci3)

        logsheet.write(11, 8, "Number of Reads (Combined)", bold)
        for i in range(0, query_length):
            logsheet.write(i + 12, 8, out_tally_all[i])
        logsheet.write(11, 9, "% of Total Reads (Combined)", bold)
        for i in range(0, query_length):
            logsheet.write(i + 12, 9, out_tally_all[i] / total, percentage_format)
        logsheet.write(11, 10, "Normalized Read Count (Combined)", bold)
        if max(out_tally_all) > 0:
            for i in range(0, query_length):
                logsheet.write(i + 12, 10, out_tally_all[i] / max(out_tally_all), deci3_format)

        logsheet.write(11, 11, "Example Reads RL", bold)
        for i in range(0, query_length):
            logsheet.write(i + 12, 11, str(example_reads_rl[i]), red_font)

        logsheet.write(11, 12, "Example Reads LR", bold)
        for i in range(0, query_length):
            logsheet.write(i + 12, 12, str(example_reads_lr[i]), blue_font)

        # 'highlight box', take from top_3 list determined above
        logsheet.write(0, 4, 'Most Frequent Transposition Distances (bp)', bold)

        logsheet.write(1, 4, top_3[0][1])
        logsheet.write(1, 5, top_3[0][0] / total, percentage_format)
        logsheet.write(1, 6, out_tally_rl[top_3[0][1]]/total, red_percent)
        logsheet.write(1, 7, out_tally_lr[top_3[0][1]] / total, blue_percent)

        logsheet.write(2, 4, top_3[1][1])
        logsheet.write(2, 5, top_3[1][0] / total, percentage_format)
        logsheet.write(2, 6, out_tally_rl[top_3[1][1]] / total, red_percent)
        logsheet.write(2, 7, out_tally_lr[top_3[1][1]] / total, blue_percent)

        logsheet.write(3, 4, top_3[2][1])
        logsheet.write(3, 5, top_3[2][0] / total, percentage_format)
        logsheet.write(3, 6, out_tally_rl[top_3[2][1]] / total, red_percent)
        logsheet.write(3, 7, out_tally_lr[top_3[2][1]] / total, blue_percent)

        logsheet.write(4, 4, 'On Target Frequency', bold)
        logsheet.write(4, 5, on_target_total / total, percentage_format)

        logsheet.write(5, 4, 'Orientation Bias (R->L:L->R)', bold)
        logsheet.write(5, 5, '{} : 1'.format(round(len(final_list_rl)/(len(final_list_lr)+0.00000001), 2)))  # in case LR is 0

    # plot and save graphs if plot setting = True
    # Only plots a certain window (e.g. from 40bp to 60bp) - change that in the xlim options below
    if plot:
        x_axis = []  # artificial x-axis
        for i in range(20, 61):
            x_axis.append(i)
        y_rl = out_tally_rl[20:61]
        y_lr = out_tally_lr[20:61]

        max_y = max(max(y_rl), max(y_lr))  # for scaling y axis
        if not plot_overlap:
            fig, axs = plt.subplots(1, 2)
            fig. tight_layout(rect=[0.15, 0.1, 1, 0.9])
            title = fig.suptitle("{} - {}\nOn-target frequency: {}%\nOrientation bias (R->L:L->R): {}:1"
                                 .format(code, description, round(100*on_target_total/total, 1),
                                 round(len(final_list_rl)/(len(final_list_lr)+0.00000001), 2)))
            title.set_y(0.88)
            # first graph on the left in red
            axs[0].bar(x_axis, y_rl, color='tab:orange', width=1.0)
            axs[0].set_title("R->L Integration Events")
            # second graph on the right in blue
            axs[1].bar(x_axis, y_lr, color='tab:blue', width=1.0)
            axs[1].set_title("L->R Integration Events")
            fig.subplots_adjust(wspace=0.7)

            for axs in axs.flat:
                axs.spines['top'].set_visible(False)
                axs.spines['right'].set_visible(False)
                # ax.spines['bottom'].set_visible(False)
                # axs.spines['left'].set_visible(False)
                axs.spines['bottom'].set_position('zero')
                axs.spines['left'].set_bounds(0, max_y)
                axs.set_xticks([40,45,50,55,60])
                axs.set_xticklabels([40,45,50,55,60])
                axs.set_yticks([0, max_y])
                axs.set_yticklabels([0, max_y])
                axs.set_xlim(left=39, right=61) ## Change window here
                axs.set_ylim(bottom=0, top=1.05*(max_y))
                axs.set(xlabel="Distance from target site (bp)", ylabel="Read count")
                axs.yaxis.set_label_coords(-0.1,0.5)

            fig.set_size_inches(6, 4.2)
            fig.subplots_adjust(top=0.65)
            #plt.xlabel("Distance from target site (bp)")
            #plt.ylabel("Read count")
            #plt.gca().set_xlim(left=35, right=60)
            #plt.gca().set_ylim(bottom=0, top=(1.3*max_y))

            # plt.savefig('test.svg', dpi=250)
            plt.savefig('Dist_Output\\{}_{}_{}.png'.format(date, code, description), dpi=300)
            plt.close()
        if plot_overlap:
            fig, axs = plt.subplots(1, 1, tight_layout=True)
            title = fig.suptitle("{} - {} / On-target = {}% / Bias = {} :1".format(
                code, description, round(100*on_target_total/total, 1),
                                 round(len(final_list_rl)/(len(final_list_lr)+0.00000001), 2)))
            title.set_y(0.9)
            # LR graph is colorless with a border
            axs.bar(x_axis, y_lr, color='none', edgecolor='#153C6B', linewidth=1.0, width=1.01, zorder=1)
            # RL graph is blue with no border (behind bordered RL)
            axs.bar(x_axis, y_rl, color='#83B0DD', edgecolor='#83B0DD', linewidth=1.0, width=1.01, zorder=0)

            axs.spines['top'].set_visible(False)
            axs.spines['right'].set_visible(False)
            # ax.spines['bottom'].set_visible(False)
            # axs.spines['left'].set_visible(False)
            axs.spines['bottom'].set_position('zero')
            axs.spines['left'].set_bounds(0, max_y)
            axs.set_xticks([40, 42, 45, 50, 55, 60])
            axs.set_xticklabels([40, 0, 45, 50, 55, 60])
            axs.set_yticks([0, max_y])
            axs.set_yticklabels([0, max_y])
            axs.set_xlim(left=42, right=58) ## Change window here
            axs.set_ylim(bottom=0, top=1.25 * (max_y))
            axs.set(xlabel="Distance from target site (bp)", ylabel="Read count")
            axs.yaxis.set_label_coords(-0.05, 0.4)

            fig.set_size_inches(5, 4.2)
            # plt.xlabel("Distance from target site (bp)")
            # plt.ylabel("Read count")
            # plt.gca().set_xlim(left=35, right=60)
            # plt.gca().set_ylim(bottom=0, top=(1.3*max_y))

            # plt.savefig('test.svg', dpi=250)
            plt.savefig('Dist_Output_Overlap\\{}_{}_overlapped_{}.svg'.format(date, code, description), dpi=500)
            plt.close()


# analyze all input files back-to-back using a master input .csv file
with open('input_spacers.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader, None)
    for row in reader:
        code = row[0]
        filename = 'none'
        for i in os.listdir('.'):
            if fnmatch.fnmatch(i, "*{}*.fastq".format(code)):
                filename = i
                break
        if filename != 'none':
            description = row[6]
            psl = row[7]
            direction = str(row[4]).lower()  ## spacer direction
            ## set refeq based on direction
            if direction == 'fw':
                refseq = genome
            elif direction == 'rv':
                refseq = genome_rc
            else:
                print("Direction Error = {}".format(code))
                break
            spacer = str(row[10]).upper()  #spacer sequence, convert to uppercase
            # run main dist_find function with the relevant input variables
            dist_find(code, psl, description, filename, direction, refseq, spacer, exp_date, excel, plot, plot_overlap)
        else:
            print("WARNING - File Not Found For {}".format(code))

if excel:
    log.close()