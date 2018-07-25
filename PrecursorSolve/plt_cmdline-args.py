## Generic Plotting Outline for rhovsenergy.out data
##
## Author:  Zander Mausolff
## Usage:  python plt_cmdline-args.py
## Usage help: python plt_cmdline-args.py help
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import re

##-------Change font to 'Palatino' throughout the plot---##
##-------This will match latex documents exactly---------##
from matplotlib import rcParams
#rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = ['Palatino']

current_dir = os.getcwd()
##-------Check if user wants help w/cmd line options
arglist = str(sys.argv)
if "help" in arglist:
	print 'user wants help'
	print 'cmd line options are: \n f - will read files for plotting from a file \n read - type file name containing data file names \n r - plot reactivity \n p - power \n np - normalize power \n  y - plot yield \n xl - x log scale \n yl - y log scale \n pdf - save figure as pdf \n png - save figure as png \n s - plot figure in window'
	exit()

twiglall= np.loadtxt(current_dir + "/" + 'twiglall_ramp.txt',skiprows=2)

f=plt.figure()

#-------Have x log scale
if re.search(r'\bxl\b',arglist[1:]):
        plt.xscale('log')
#-------Have y log scale
if re.search(r'\byl\b',arglist[1:]):
        plt.yscale('log')

if re.search(r'\bread\b',arglist[1:]):
data_file_names="data_file_names.txt"
file_names = [line.rstrip('\n') for line in open(data_file_names)]

i=0
while i < len(file_names):
	data1 = np.loadtxt(current_dir + "/" + file_names[i] , skiprows=2)
	x_coord = data1[:,0]
    prec_conc = data1[:,5]
    plt.plot(x_coord,prec_conc,label='Test')
    plt.ylabel('Precursor Concentration')
	plt.xlabel('Distance')
    i=i+1

##------Evaluate command line
length = len(str(arglist[1]))
print 'Argument List:', arglist[length:]


##------Plot benchmark data

##------Configure the legend --- ##
plt.legend(loc='lower right',prop={'size':14},numpoints=1)

if re.search(r'\bpng\b',arglist[1:]):
	f.savefig("./" + name + ".png")

if re.search(r'\bpdf\b',arglist[1:]):
        f.savefig("./" + name + ".pdf")

if re.search(r'\bs\b',arglist[1:]):
	plt.show()
