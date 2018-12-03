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

##-------Useful for animation
from matplotlib import animation
from matplotlib.animation import FuncAnimation

## Get current command line location
current_dir = os.getcwd()

##-------Check if user wants help w/cmd line options
arglist = str(sys.argv)

##------Shape functions------##
def shape_fcn(fcn,x,i):
    if (i == 0) :
        fcn = -0.5*x*(1.0 - x) 
    if (i == 1) :
        fcn = (1 + x)*(1 - x)
    if (i == 2) :
        fcn = 0.5*x*(1 + x)
    return fcn

# Configure figure
fig, ax = plt.subplots()
#fig.set_tight_layout(True)

# Read files in
data_file_names="data_file_names.txt"
file_names = [line.rstrip('\n') for line in open(data_file_names)]

#i=0
color = ['-','-','-','-']

#set number of precursor groups
num_precursor_groups = 1 

data_initial = np.loadtxt(current_dir + "/" + file_names[0] , skiprows=2)
x_coord_initial = data_initial[:,0]
max_length = len(x_coord_initial)
element_length = x_coord_initial[2] - x_coord_initial[0]
starting = x_coord_initial[0]
ending =   x_coord_initial[max_length-1]
elem_interval = 100
norm_space = np.linspace(-1,1,elem_interval)
x_interval = (ending*elem_interval)/element_length
# Break the x coordinate into smaller intervals
x_space_initial = np.linspace(starting,ending,x_interval)
line_mod, = ax.plot([],[],'-')

time_int = 0.0
x_max = 550 
y_max = 5.5E9 
def update(i):
    data1 = np.loadtxt(current_dir + "/" + file_names[i] , skiprows=2)
    g=0
    while g < num_precursor_groups:
        plt.clf()
        time_interval=0.1
        x_coord = data1[:,0]
        prec_conc = data1[:,g+3]
        max_length = len(x_coord)
        element_length = x_coord[2] - x_coord[0]
        starting = x_coord[0]
        ending =   x_coord[max_length-1]
        #elem_interval = 100
        
        norm_space = np.linspace(-1,1,elem_interval)
        x_interval = (ending*elem_interval)/element_length
        # Break the x coordinate into smaller intervals
        x_space = np.linspace(starting,ending,x_interval)
        
        # evaluate using quadratic interpolation functions
        nodes_per_elem = 3
        num_elem = max_length/nodes_per_elem
        value = 0
        
        # array to hold quadratically interpolated solution
        soln = []
        # evaluate solution between nodal points
        for q in range(0, num_elem):
            ii = (q)*nodes_per_elem 
            for k in range(0,len(norm_space)):
                x = norm_space[k]
                linear_comb = 0
                for j in range(0,nodes_per_elem):
                    value = shape_fcn(value,x,j)
                    linear_comb = linear_comb + value*prec_conc[j + ii]
                soln.append(linear_comb)  
        print file_names[i] 
        plt.ylim(0,y_max)
        plt.xlim(0,x_max)

        time_int = i*time_interval
        plt.plot(x_space_initial,soln,'--', label='time = ' + ' %.2f' % time_int + ' s')
        #line_mod.set_ydata(soln) 
        #line_mod.set_xdata(x_space_initial)
        
        plt.ylabel('Precursor Concentration',size=14)
        plt.xlabel('Distance [cm]',size=14)
        plt.legend()
        #plt.legend('time = '+ (file_names[i][len(file_names)-3:len(file_names)])) 
        g=g+1
    return line_mod, ax    

anim = FuncAnimation(fig, update, frames=np.arange(0,len(file_names)), interval=200)

#anim.save('line.gif', dpi=300, writer='imagemagick')

plt.show()

## Name of saved image
name = 'precursor_concentration'

##------Configure the legend --- ##
#plt.legend( loc='upper center', bbox_to_anchor=(0.5, -0.15),  shadow=True,prop={'size':14},numpoints=1)

if re.search(r'\bpng\b',arglist[1:]):
	f.savefig("./" + name + ".png")

if re.search(r'\bpdf\b',arglist[1:]):
        f.savefig("./" + name + ".pdf", bbox_inches = "tight")

if re.search(r'\bs\b',arglist[1:]):
	plt.show()
