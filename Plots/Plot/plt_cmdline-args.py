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
from matplotlib.font_manager import FontProperties

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

#twiglall= np.loadtxt(current_dir + "/" + 'twiglall_ramp.txt',skiprows=2)

f=plt.figure()

#-------Have x log scale
if re.search(r'\bxl\b',arglist[1:]):
        plt.xscale('log')
#-------Have y log scale
if re.search(r'\byl\b',arglist[1:]):
        plt.yscale('log')

if re.search(r'\bf\b',arglist[1:]):
	if re.search(r'\bread\b',arglist[1:]):
		data_file_names = raw_input("Please enter the name of the file containg the data file names:")
		print "Input names are contained in file:", data_file_names
	else:
		data_file_names="data_file_names.txt"
	file_names = [line.rstrip('\n') for line in open(data_file_names)]

	i=0
	while i < len(file_names):
		data1 = np.loadtxt(current_dir + "/" + file_names[i] , skiprows=1)
		time1 = data1[:,0]
		if re.search(r'\br\b',arglist[1:]):
			reactivity1 = data1[:,5]
			plt.plot(time1,reactivity1,label='RK4 Solver')
			plt.ylabel('Reactivity [$]')
			plt.xlabel('Time [s]')
			name="reactivity"
		if re.search(r'\bnp\b',arglist[1:]):
			plt.ylabel('Relative Power',fontsize=14)
			plt.xlabel('Time [s]',fontsize=14)
			name="norm_power"
			power1 = data1[:,3]
			pow_first=power1[0]
			norm_power=power1/pow_first
			plt.plot(time1,norm_power,label='RK4 solver')
		if re.search(r'\bp\b',arglist[1:]):
			plt.ylabel('Power ')
			plt.xlabel('Time [s]')
			name="power"
			power1 = data1[:,1]
			pow_first=power1[0]
			norm_power=power1/pow_first
			plt.plot(time1,power1,label=file_names[i][25:len(file_names[i])])
		if re.search(r'\by\b',arglist[1:]):
			plt.ylabel('Precursor Concentration')
			plt.xlabel('Time [s]',fontsize=14)

			name="precursor_conc"
			precursor1 = data1[:,2]
			precursor2 = data1[:,3]
			precursor3 = data1[:,4]
			precursor4 = data1[:,5]
			precursor5 = data1[:,6]
			precursor6 = data1[:,7]
			precursor7 = data1[:,8]
			precursor8 = data1[:,9]

			plt.plot(time1,precursor1,label='precursor 1')
			plt.plot(time1,precursor2,label='precursor 2')
			plt.plot(time1,precursor3,label='precursor 3')
			plt.plot(time1,precursor4,label='precursor 4')
			plt.plot(time1,precursor5,label='precursor 5')
			plt.plot(time1,precursor6,label='precursor 6')
			plt.plot(time1,precursor7,label='precursor 7')
			plt.plot(time1,precursor8,label='precursor 8')

		i=i+1
else:
	print "You are not reading file names from a text file, you must manually edit this python file"
##------Evaluate command line

length = len(str(arglist[1]))
print 'Argument List:', arglist[length:]

fontP = FontProperties()
fontP.set_size('small')


##------Plot benchmark data

##------Configure the legend --- ##
#plt.legend(loc='lower right',prop={'size':14},numpoints=1)
lgd = plt.legend(loc='upper center',bbox_to_anchor=(0.5, -0.2),  shadow=True, ncol=2)

if re.search(r'\bpng\b',arglist[1:]):
	f.savefig("./" + name + ".png",bbox_extra_artists=(lgd,), bbox_inches='tight')

if re.search(r'\bpdf\b',arglist[1:]):
        f.savefig("./" + name + ".pdf",bbox_extra_artists=(lgd,), bbox_inches='tight')

if re.search(r'\bs\b',arglist[1:]):
	plt.show()
