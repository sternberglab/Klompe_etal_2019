py_name = 'bin_reindex_v1.py'

#Reindexes the geneious coverage csv so that the end of the protospacer is at the start of its respective bin
#Returns a new reindexed coverage file

import os
import csv
import fnmatch

os.chdir('C:\\Users\\Leo Vo\\Desktop\\SEK_NGS\\NGS data cleaned unique hits\\Vch - Fastq\\Trans_Sites')

bin_size = 100

mode = 'Trans_Sites'   # 'Trans_Sites' or 'Coverage'
# However this does not change the end of the name of the output file, just be careful.
with open('indv_test.csv', newline='') as csvfile:  #open master csv file containing names of files to reindex
    reader = csv.reader(csvfile, delimiter=',')
    next(reader, None)  # skips one line (header)
    for row in reader:
        code = row[0]
        filename = 'none'
        for i in os.listdir('.'):  # look for files matching name
            if fnmatch.fnmatch(i, "{}*{}.csv".format(code, mode)):
                filename = i
                break
        if filename != 'none': # proceed if file is found
            spacer = int(row[8]) # get spacer position from master csv
            start = 0 # default start is 0
            if str(row[4]).lower() == 'fw':
                start = spacer - (int(spacer/bin_size)*bin_size) # change start to spacer
            elif str(row[4]).lower() == 'rv':
                start = spacer - (int(spacer/bin_size)*bin_size) # change start to spacer
            else:
                print("WARNING - No Direction")
            with open(filename, newline='') as csv_in: # open csv file containing original coverage values
                reader_in = csv.reader(csv_in, delimiter=',') # read out these values
                in_lines = list(reader_in) # put them in a list in order
                in_lines.remove(in_lines[0]) # remove first value (header)
            csv_in.close()
            with open('REbin_{}_{}bp_Coverage.csv'.format(code, bin_size), 'w', newline='') as csv_out: # open new csv
                writer = csv.writer(csv_out)
                writer.writerow(['Position_rebinned', 'Coverage']) # make header line
                for i in range(start, len(in_lines)): # begin writing from the new start position, to the end of list
                    writer.writerow(in_lines[i])
                for i in range(0, start): # then write in values from beginning of list, to the new start position
                    writer.writerow(in_lines[i])
            csv_out.close()

