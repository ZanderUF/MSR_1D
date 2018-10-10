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
from matplotlib import pylab

##-------Change font to 'Palatino' throughout the plot---##
##-------This will match latex documents exactly---------##
from matplotlib import rcParams

current_dir = os.getcwd()
##-------Check if user wants help w/cmd line options
arglist = str(sys.argv)

f=plt.figure()

##------Shape functions------##
def shape_fcn(fcn,x,i):
    if (i == 1) :
        fcn = -0.5*x*(1.0 - x) 
    if (i == 2) :
        fcn = (1 + x)*(1 - x)
    if (i == 3) :
        fcn = 0.5*x*(1 + x)
    return fcn 

# Create space to evaluate shape functions
elem_interval = 100
norm_space = np.linspace(-1,1,elem_interval)

# evaluate using quadratic interpolation functions
nodes_per_elem = 3
value = 0

shape_1=0
shape_2=0
shape_3=0
# array to hold quadratically interpolated solution
shape_array_1 = []
shape_array_2 = []
shape_array_3 = []
# evaluate between nodal points
for k in range(0,len(norm_space)):
    x = norm_space[k]
    shape_1 = shape_fcn(shape_1,x,1)
    shape_2 = shape_fcn(shape_2,x,2)
    shape_3 = shape_fcn(shape_3,x,3)
    shape_array_1.append(shape_1)  
    shape_array_2.append(shape_2)
    shape_array_3.append(shape_3)

plt.plot(norm_space,shape_array_1,'r-', label='$f_e^1(x)$')
plt.plot(norm_space,shape_array_2,'b.', label='$f_e^2(x)$')
plt.plot(norm_space,shape_array_3,'g--', label='$f_e^3(x)$')
plt.axhline(0,linewidth=0.8, color='black')

plt.ylabel('$f_e(x)$',size=14)
plt.xlabel('x', labelpad = 1, size=14)

## Name of saved image
name = 'solution_approx'

##------Configure the legend --- ##
pylab.legend(loc=9, bbox_to_anchor=(0.5, -0.15), ncol=3)

if re.search(r'\bpng\b',arglist[1:]):
	f.savefig("./" + name + ".png",bbox_inches="tight")

if re.search(r'\bpdf\b',arglist[1:]):
        f.savefig("./" + name + ".pdf",bbox_inches="tight")

if re.search(r'\bs\b',arglist[1:]):
	plt.show()
