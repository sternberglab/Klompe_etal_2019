py_name = 'linear_hist_sites_batch_v1B.py'

################
# Leo Vo Feb 2019

# Script Function:
# Bins coverage/integration read data across the genome additively into uniform bins of defined size. Outputs a bar
# plot for visualization.
# v1B: input .csv files now count duplicated reads

# Script Input(s):
# 1/ Csv files from either Geneious output (as Coverage) or from python code (as Reads representing
# integration 'events'), containing coordinates and count (2 columns). These csvs needs to have exactly 1 header row.

# 2/ One csv file containing information about the csv files being run (such as Sample ID, Description, etc)
# In order to run multiple csv files in one batch continuously.

# Script Output(s):
# One genome-wide Tn-Seq bar plot with peaks corresponding to where Coverage/Integration Reads are found, per input csv.
#################


import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import fnmatch

# font control
plt.rcParams['svg.fonttype'] = 'none'  # important so that text stays as intact characters in the output
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

# change directory based on where input files are
os.chdir('C:\\Users\\Leo Vo\\Desktop\\SEK_NGS\\190217_Integration_site_analysis_SEK\\Trans_Sites')

# turns of interactive mode for matplotlib, optional
plt.ioff()

# some run parameters
bin_size = 5000  # in bp
exp_date = 190202
mode = 'Trans_Sites'  # either 'Coverage' or 'Trans_Sites', for csv input file discovery
norm_on = True  # normalize each data point against total (genome mapped) reads in each sample (as percentage)

def readnbin(filename, bin_size):  # main binning function, edited for duplicated reads counting
    # read in csv file, extract total from first line
    with open(filename, newline='') as f:
        reader = csv.reader(f)
        row1 = next(reader)
        total = int(row1[1])  # get total reads, not the same as sum of all coverages (due to duplication)

    # extract the rest of the coverage values
    csvinput = np.genfromtxt(filename, delimiter=",", skip_header=2)  # skips 2 header rows
    z = csvinput[:,1]  # takes only the "Coverage" column into the input list
    a = [int(i) for i in z]  # convert everything from floats to int

    # sets up empty list which will become the output list
    a2_size = int((len(a)-1)/bin_size) + 1  # determines size of out list (number of bins post binning)
    a2 = [0]*a2_size

    # bin data additively
    y = 0
    while y <= len(a):
        x = 0
        while x < bin_size:
            if y + x < len(a):
                a2[int(y / bin_size)] += a[y + x]
                x += 1
            else:
                x += bin_size
        y += bin_size

    # set up second list that will serve as x-axis ('coordinate') for output graph
    b = []
    for i in range(1, len(a2)+1):
        b.append(i)

    # convert lists to arrays for matplotlib
    b = np.asarray(b)
    a2 = np.asarray(a2)

    return b, a2, total  # return the lists. If binning using coverage, also return max(a) here for ytick labels later


def plot(filename, bin_size, date, desc, psl, spacer_location, norm_on):  # main graphing function
    spacer_bin = int(spacer_location/bin_size) + 0.5  # determine which bin the spacer lies in

    # set up figure and call readnbin()
    fig, axs = plt.subplots(1, 1, tight_layout=True)
    b, a2, total_reads = readnbin(filename, bin_size)  # call readnbin to get arrays for xy axes

    if norm_on:  # normalize based on total reads (percentage)
        a2 = (100*a2)/total_reads  # make sure a2 is returned as an array from readnbin()

    # change y-max based on normalization option
    if norm_on:
        max_y = 100
    else:
        max_y = max(a2)  # if using coverage, take max(a) from readnbin() instead of max(a2) (the 'pre binning' maximum)

    axs.scatter(spacer_bin, -0.03 * max_y, marker="^", c="#A32E79", s=33)  # 1 scatter plot point for protospacer marker
    axs.bar(b, a2, color='#0D62AC', width=1.05)  # main bar graph, width slightly >1 to avoid gaps between bars

    # configure spines of the graph; 'remove' top and right spines
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.spines['bottom'].set_position('zero')
    axs.spines['left'].set_bounds(0, max_y)
    axs.spines['bottom'].set_bounds(0, len(a2))  # bottom line starts and end with length of genome
    axs.spines['left'].set_position('zero')

    # 2 ticks on the y axis, one at 0 and on at max value of y (max reads)
    axs.set_yticks([0, max_y])

    if norm_on:
        axs.set_yticklabels([0, '100%'])
    else:
        axs.set_yticklabels([0, max(a2)])

    # set up xticks
    bin_scale = int(100000/bin_size)
    axs.set_xticks([0,5*bin_scale,10*bin_scale,15*bin_scale,20*bin_scale,25*bin_scale,30*bin_scale,35*bin_scale,
                    40*bin_scale,45*bin_scale])
    axs.set_xticklabels([0,0.5,1,1.5,2,2.5,3,3.5,4,4.5])

    # size of output figures
    fig.set_size_inches(8.9, 3.1)

    # axis labeling
    #plt.xlabel("E. coli genomic coordinate (Mbp)")
    #plt.ylabel("Read count")
    axs.yaxis.set_label_coords(-0.02, 0.5)

    # set limits of graph area - different from spine limits
    # leave extra space on top for info text, and space on bottom for the spacer marker
    plt.gca().set_xlim(left=-0.15*len(a2), right=len(a2))
    plt.gca().set_ylim(bottom=-0.12*max_y, top=1.20*max_y)

    # info text
    plt.text(5*bin_scale, 1.17*max_y, "{}-{}\n{}\ngRNA Plasmid = {}\nTotal Reads = {}"
             .format(code, date, desc, psl, total_reads))

    # save figures
    plt.savefig('plots\\DUP_{}_{}_{}.svg'.format(date, code, desc), dpi=500)
    #plt.savefig('plots\\{}_{}_{}.png'.format(date, code, desc), dpi=500)
    plt.close()  # closes the matplotlib preview popup window


# read in the info csv file and begin processing files as a batch
with open('input_spacers_test.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader, None)  # skips header row
    for row in reader:
        code = row[0]
        # search for file in cwd based on Sample ID ('code')
        filename = 'none'
        for i in os.listdir('.'):
            if fnmatch.fnmatch(i, "{}*{}.csv".format(code, mode)):
                filename = i
                break
        if filename != 'none':
            desc = row[6]
            psl = row[7]
            if row[8] == '-' or row[8] == ' ':
                spacer_location = -100000000  # set a large negative value if spacer isn't found
            else:
                spacer_location = int(row[8])
            plot(filename, bin_size, exp_date, desc, psl, spacer_location, norm_on)
        else:
            print("WARNING - File Not Found For {}".format(code))

