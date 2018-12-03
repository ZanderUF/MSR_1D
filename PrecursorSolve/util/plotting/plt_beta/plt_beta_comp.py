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

current_dir = os.getcwd()
##-------Check if user wants help w/cmd line options
arglist = str(sys.argv)

f=plt.figure()

##------Shape functions------##
def shape_fcn(fcn,x,i):
    if (i == 0) :
        fcn = -0.5*x*(1.0 - x) 
    if (i == 1) :
        fcn = (1 + x)*(1 - x)
    if (i == 2) :
        fcn = 0.5*x*(1 + x)
    return fcn 

data_file_names="data_file_names.txt"
file_names = [line.rstrip('\n') for line in open(data_file_names)]
tick_marks = ['--',':','-','.']

i=0
color = ['.', 'o', '*', '^']

#set number of precursor groups
labels = ['Lagged','Instant']
skp=1
while i < len(file_names):
    print file_names[i]

    data1 = np.loadtxt(current_dir + "/" + file_names[i] , skiprows=skp)
    time = data1[:,0]
    mass_flow = data1[:,5]
    beta   = data1[:,4]

    plt.plot(time,beta,tick_marks[i],label=labels[i] )
    plt.ylabel('Beta',size=14)
    
    #plt.xlabel('Mass Flow [g/s]',size=14)
    plt.xlabel('Time [s]',size=14)
    plt.ylim(6.00E-3, 7.0E-3)    
    i=i+1

## Name of saved image
name = 'beta_variation'

##------Configure the legend --- ##
plt.legend( loc='upper left',  shadow=True)

if re.search(r'\bpng\b',arglist[1:]):
	f.savefig("./" + name + ".png", bbox_inches = "tight")

if re.search(r'\bpdf\b',arglist[1:]):
        f.savefig("./" + name + ".pdf", bbox_inches = "tight")

if re.search(r'\bs\b',arglist[1:]):
	plt.show()
