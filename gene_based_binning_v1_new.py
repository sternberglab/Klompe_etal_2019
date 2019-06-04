py_name = 'gene_based_binning_v1_new.py'

################
# Leo Vo Feb 2019

# Script Function:
# Bin coverage of reads into coding sequence ('gene') intervals. Then plot 2d replicate analysis graph
# Mostly used for analysis of random transposition (e.g. Himar)

# Some details:
# define bins based on intervals (from ecoli_genes.txt, which was made with Eco_Seq_Gene_Parse.py):
    ## extract the genes file into a list of lines; intervals are lines[i].split(' ')[2] to [3]
# bin at position 0 in the tally list is interval 1 (from line 0 in the file)
# number of bins is the number of lines in the .txt input, or number of genes
# take coordinates and check to see whether it lies between intervals, then increase counter in tally list.

# Script Input(s):
# 1/ ecoli_genes_txt, containing coordinates (start & end) of coding sequences
# 2/ .csv files containing tallied transposition coordinates of reads

# Script Output(s):
# 1/ A 2D plot of the value of each bin within each of 2 replicates
# 2/ Optional: .csv outputs of binned reads (single column) for separate statistical analyses
#################

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv
import matplotlib.ticker as ticker

# change directory based on where input files are
os.chdir('C:\\Users\\Leo Vo\\Desktop\\SEK_NGS\\NGS data cleaned unique hits\\Himar\\Trans_Sites')
# define input .csv files, in order of rep 1, rep 2, etc.
files = ['A4796_Trans_Sites.csv', 'A4788_Trans_Sites.csv', 'A4781_Trans_Sites.csv', 'A4754_Trans_Sites.csv']

# font control
plt.rcParams['svg.fonttype'] = 'none'  # important so that text stays as intact characters in the output
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

# some run parameters
log_on = True  # plot on a log scale (also affects csv output)
log_scale = 2
save_csv = True  # output csvs of binned data for statistical analyses
lin_reg = True  # draw linear regression line on graph (have to input slope etc. manually)

# read in parsed gene intervals file and put all lines into a list
with open('ecoli_genes.txt', 'r') as open_file:
    lines = open_file.readlines()


def gene_bin(filename, lines):  # main binning function
    # read in csv file and take "read count/coverage/events" column
    csv = np.genfromtxt(filename, delimiter=",", skip_header=1)
    z = csv[:, 1]
    a = [int(i) for i in z]  # convert everything into integers

    # set up tally list based on number of genes
    tally = [0]*len(lines)

    # each coding sequence interval corresponds to a bin in tally
    # split() function returns the 2 coding sequence coordinates
    # loop through a in between each interval, sum every 'coverage' number into the corresponding bin in tally

    for gene in range(0, len(lines)):
        for i in range(int(lines[gene].split(' ')[2]), int(lines[gene].split(' ')[3]) + 1):  # +1 since inclusive
                tally[gene] += a[i]
    return tally

# call main function for each rep
# rep1 = gene_bin(files[0], lines, log_on, log_scale)
rep2 = gene_bin(files[1], lines)
rep3 = gene_bin(files[2], lines)

# Compare the reps together and remove any 0 points
def remove_zero_log(a, b, log_scale):  # remove (0,0) points, use for log scale
    # Use -1 so that there is a new position for the (positive,0) points.
    # Will need to relabel the axis markers to 0, 2^0, 2^1 etc..
    a_new = []
    b_new = []
    for i in range(0, len(a)):
        if a[i] > 0 and b[i] > 0:
            a_new.append(math.log(a[i], log_scale))
            b_new.append(math.log(b[i], log_scale))
    return a_new, b_new

if log_on:
    rep23, rep32 = remove_zero_log(rep2, rep3, log_scale)


# brief function for csv output
def out_csv(file, data):
    with open('test_new_{}_Binned.csv'.format(file), 'w', newline='') as csv_out:
        writer = csv.writer(csv_out)
        writer.writerow(['Reads'])  # include header row
        for i in range(0, len(data)):
            writer.writerow([data[i]])
    csv_out.close()


# output CSVs for stats analyses:
if save_csv:
    if log_on:
        out_csv('A4788', rep23)
        out_csv('A4781', rep32)
    else:
        out_csv('A4788', rep2)
        out_csv('A4781', rep3)
# initialize figure
fig = plt.figure()

# change based on which reps are used
if log_on:
    x = rep23
    y = rep32
else:
    x = rep2
    y = rep3

axs = fig.add_subplot(111)
axs.scatter(x, y, s=0.5, c='k', marker='o')
axs.set_title("Himar, Rep-1-2, gene-based binning")

# spines and labels
axs.spines['top'].set_visible(False)
axs.spines['right'].set_visible(False)
axs.set_xlabel('Log2 Reads, Replicate 1')
axs.set_ylabel('Log2 Reads, Replicate 2')
plt.gca().set_xlim(left=-0.1*max(x), right=1.1*max(x))
plt.gca().set_ylim(bottom=-0.1*max(y), top=1.1*max(y))

# control tick frequency
#tick_spacing = 2
#axs.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
#axs.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

if lin_reg:  # plots linear regression line based on manually inputted parameters
    intercept = 1.0827857234699376
    slope = 0.8647442059333533
    x_vals = np.array(axs.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals)

plt.tight_layout()
fig.set_size_inches(6,6)
plt.savefig('2d\\retry_new_Himar_23_linear_genebin_log2_linregress.svg', dpi=500)  # make sure subfolder '2d\' exists
#plt.show()
plt.close()  # disable pop up graph output

###
count3 = 0
for i in rep3:
    if i == 0 or i == 1:
        count3 += 1

count32 = 0
for i in rep32:
    if i == 0:
        count32 += 1