py_name = 'coverage_rep_analysis_v1_new.py'

# Leo Vo Feb 2019

# Script Function:
# Plots a 2D graph of integration reads/events to probe for reproducibility of genome wide read coverage across

# Script Input(s):
# 1/ Csv files from either Geneious output (as Coverage) or from python code (as Reads representing
# integration 'events'), containing coordinates and count (2 columns). These csvs needs to have exactly 1 header row.


# Script Output(s):
# 1/ A 2D plot of the value of each bin within each of 2 replicates. If plotting in log scale, will need to manually
# relabel the axis markers because of how (positive, 0) points are plotted.


from mpl_toolkits.mplot3d import Axes3D  # have here if doing 3D plots with matplotlib
import os
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.ticker as ticker

# change directory based on where input files are
os.chdir('C:\\Users\\Leo Vo\\Desktop\\SEK_NGS\\NGS data cleaned unique hits\\Vch - Fastq\\Trans_Sites')

# font control
plt.rcParams['svg.fonttype'] = 'none'  # important so that text stays as intact characters in the output
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

# run parameters and define input files (as rep 1, rep 2, etc.)
bin_size = 100  # size of bins in bp
log_on = True  # to turn on log scale
log_scale = 2

#dump input files in the list below, call them later by using list indices
#files = ['REbin_A4738_{}bp_Coverage.csv'.format(bin_size),
#         'REbin_A4739_{}bp_Coverage.csv'.format(bin_size),
#         'REbin_A4756_{}bp_Coverage.csv'.format(bin_size)]
#files = ['A4783_Trans_Sites.csv', 'A4749_Trans_Sites.csv', 'A4759_Trans_Sites.csv']
files = ['A4783_Trans_Sites.csv', 'A4748_Trans_Sites.csv', 'A4750_Trans_Sites.csv']

def readnbin(filename, bin_size):  # main binning function
    # read in csv files
    csv = np.genfromtxt(filename, delimiter=",", skip_header=1)
    z = csv[:,1]  # takes only the "Coverage" column into the input list
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

    # convert list to arrays for matplotlib, optional??
    a2 = np.asarray(a2)

    return a2  # return a2. (only 1 list needed here). One input csv returns one list


def remove_zero(a, b):  # remove (0,0) points, use for linear scale
    a_new = [0]*len(a)
    b_new = [0]*len(b)
    # next three lists are used for labeling points that are >0 for both axes later
    noz_bin_nums = []
    a_noz_list = []
    b_noz_list = []
    for i in range(0, len(a)):
        if a[i] == 0 and b[i] == 0: # "remove" (0,0) points by making them very negative
            a_new[i] = -100000
            b_new[i] = -100000
        if a[i] > 0 and b[i] > 0: # add (pos,pos) points into lists above
            noz_bin_nums.append(i+1)
            a_noz_list.append(a[i])
            b_noz_list.append(b[i])
        if a[i] > 0:
            a_new[i] = a[i]
        if b[i] > 0:
            b_new[i] = b[i]
    return a_new, b_new, noz_bin_nums, a_noz_list, b_noz_list


def remove_zero_log(a, b, log_scale):  # remove (0,0) points, use for log scale
    # Use -1 so that there is a new position for the (positive,0) points.
    # Will need to relabel the axis markers to 0, 2^0, 2^1 etc..
    a_new = [-1]*len(a)
    b_new = [-1]*len(b)
    # next three lists are used for labeling points that are >0 for both axes later
    noz_bin_nums = []
    a_noz_list = []
    b_noz_list = []
    for i in range(0, len(a)):
        if a[i] == 0 and b[i] == 0:
            a_new[i] = -5000
            b_new[i] = -5000
        if a[i] > 0 and b[i] > 0:
            noz_bin_nums.append(i+1)
            a_noz_list.append(math.log(a[i], log_scale))
            b_noz_list.append(math.log(b[i], log_scale))
        if a[i] > 0:
            a_new[i] = math.log(a[i], log_scale)
        if b[i] > 0:
            b_new[i] = math.log(b[i], log_scale)
    return a_new, b_new, noz_bin_nums, a_noz_list, b_noz_list


# bin data from input csvs. **Careful with indexing here to match lists of input files above
rep1 = readnbin(files[0], bin_size)
rep2 = readnbin(files[1], bin_size)
rep3 = readnbin(files[2], bin_size)

# remove zeros from binned data (and convert to log if log_on = True)
if log_on:
    noz_12, noz_21, noz_bin_nums_12, noz_12_list, noz_21_list = remove_zero_log(rep1, rep2, log_scale)
    noz_13, noz_31, noz_bin_nums_13, noz_13_list, noz_31_list = remove_zero_log(rep1, rep3, log_scale)
    noz_23, noz_32, noz_bin_nums_23, noz_23_list, noz_32_list = remove_zero_log(rep2, rep3, log_scale)
if not log_on:
    noz_12, noz_21, noz_bin_nums_12, noz_12_list, noz_21_list = remove_zero(rep1, rep2)
    noz_13, noz_31, noz_bin_nums_13, noz_13_list, noz_31_list = remove_zero(rep1, rep3)
    noz_23, noz_32, noz_bin_nums_23, noz_23_list, noz_32_list = remove_zero(rep2, rep3)

# define xy axes based on which pair of replicates are run (noz_xy and noz_yx)
x = noz_32
y = noz_23

# set up matplotlib figure
fig = plt.figure()
axs = fig.add_subplot(111)

# optional for plotting 3D graphs:
#axs = fig.add_subplot(111, projection='3d')
#axs.scatter(x, y, z, c='r', marker='o')

axs.scatter(x, y, s=0.5, c='k')  # main scatter plot
axs.set_title("NoCascade, 100bp bins, Rep-1-3, log scale = {}".format(log_on))  # title

# axis labeling
if log_on:
    axs.set_xlabel('Reads, Replicate 1')
    axs.set_ylabel('Reads, Replicate 3')

if not log_on:
    axs.set_xlabel('Reads, Replicate 1')
    axs.set_ylabel('Reads, Replicate 3')


# make top and right spines invis
axs.spines['top'].set_visible(False)
axs.spines['right'].set_visible(False)

# (optional) label points positive for both axes with the appropriate bin numbers
for i, bin_num in enumerate(noz_bin_nums_23):
    axs.annotate(bin_num, (noz_23_list[i], noz_32_list[i]))

# graph area limits
if max(x) != 0:
    plt.gca().set_xlim(left=-1.2, right=1.1*max(x))
else:
    plt.gca().set_xlim(left=-0.1 * 1, right=1.1 * 1)
if max(y) != 0:
    plt.gca().set_ylim(bottom=-1.21, top=1.1*max(y))
else:
    plt.gca().set_ylim(bottom=-0.1 * 1, top=1.1 * 1)

# control tick frequency
tick_spacing = 1
axs.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
axs.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
#plt.gca().set_zlim(bottom=-0.1*max(z), top=1.1*max(z))

# misc plot parameters and save output
plt.tight_layout()
fig.set_size_inches(6,6)
plt.savefig('new_2d\\new_new_NoCascade_rep13_{}bp_log2.svg'.format(bin_size), dpi=500)  # make sure 2d\ subfolder exists
#plt.savefig('2d\\TnsABCQ_rep23_{}bp_log2.png'.format(bin_size), dpi=250)
plt.close()


