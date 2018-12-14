.. _Input: 

Input File Guide
================

The structure of the input files is broken up into blocks with character sequences seperating the blocks
The first line in the file is a comment, this should be something meaningful to distinguish the problem from others

.. Parameter block::

Parameter Block
---------------

Begins with **read parm** and ends with **end parm**

=========  ===========  ==========================================================
Parameter  Value        Description
=========  ===========  ==========================================================
time       yes or no    Perform a time dependent calculation
dbg        yes or no    Turn on DEBUG printing - prints a lot of stuff
meth       integer      0 - forward Euler 1 - backward Euler       
step       yes or no    Perform a step perturbation or not
stpb       Double       Start time for step perturbation
stpe       Double       End time for step perturbation
rmpb       Double       Start time for ramp perturbation
rmpe       Double       End time for ramp perturbation
ramp       yes or no    Perform a ramp perturbation or not
zag        yes or no    Perform a zig-zag perturbation or not
feed       Integer      0 - (default) no feedback 1 - temperature 2 - instant beta
rdpw       yes or no    Read power profile from a file
del        Double       Indicate fixed time step size, ex: 1E-5
tmax       Real         Indicate the total length of the simulation
tin        Double       Beginning point of the perturbation 
nem        Integer      Number of equal spaced elements in the model
npe        Integer      Number of nodes per element 
pipe       Integer      Number of elements in the outer circuit
inlt       Integer      Starting element of the inlet plenum
scor       Integer      Starting element of the main core region
ecor       Integer      Final element of the main core region
hexs       Integer      Heat exchanger starting element location
hexe       Integer      Heat exchanger ending element location
outl       Integer      Final element of the outlet plenum
area       Real         Area of the fuel core [:math:`cm^2`] 
apip       Real         Area of the piping [:math:'cm^2']
mflow      Real         Mass flow rate [:math:`g/cm^3`]
tpow       Real         Total axial power [:math:`n/s-cm`]
nitr       Real         Set number of maximum nonlinear iterations 
elem       Real         Element length [:math:`cm`]
ndg        Integer      Number of delayed neutron groups
nmat       Integer      Number of fissional materials
gen        Real         Neutron generation time [:math:`s^-2`]
reac       Real         Reactivity insertion for step and ramp
save       Real         Interval to write out spatial solution files to
=========  ===========  ==========================================================

.. Delay Block::

Delay Block
-----------

Begins with **read delay** and ends with **end delay**

=========  ===========  ================================================
Parameter  Value        Description
=========  ===========  ================================================
mat        Integer      Material identifier, not be greater than nmat
alam       Real(ndg)    Decay constant for each delayed group 
beta       Real(ndg)    Delayed neutron fraction for each delayed group
=========  ===========  ================================================

.. Sample Input File::

Sample Input File
-----------------

.. code-block:: guess

    This is a comment line
    read parm
       time=yes
       step=no
       ramp=no
       zag=yes
       del=1E-4
       tmax=10
       tin=0.0 
       nem=10
       npe=3
       pipe=6
       area=10000.0
       apip=1000.0
       mflow=5000.0
       tpow=20
       nitr=300
       elem=1.0
       ndg=6
       nmat=1
       gen=5E-4
       reac=0.0
    end parm
    read delay
      mat=1
    alam=0.0127 0.0317 0.115 0.311 1.4 3.87 end
    beta=2.85E-4 1.5975E-3 1.41E-3 3.0525E-3 9.6E-4 1.95E-4  end
    end delay

.. Input File DIF3D Values::

Input file structure to read in from DIF3D
------------------------------------------

Values from DIF3D can be read in and projected onto the domain in the 1D problem
Right now it reads in the power, power fraction, doppler reactivity worth, and density reactivity worth. These are all spatially dependent.  These values are integrated across the AREA specificed in the DIF3D problem. 

The file read in is assumed to be titled: dif3d_values.txt
It assumes the file is read in as follows

====== ================================================ =========================
Line   Value                                            Units
====== ================================================ =========================
1      Total number of axially integrated values        [integer value]
2      Spatially integrated power to be read in         [Watts]
3      Spatially integrated fraction of power           [Normalized by total] 
4      Spatially integrated doppler worth               [reactivity]
5      Total change in temperature for the perturbation [K]
6      Spatially integrated density worth               [reactivity] 
7      Total change in density during perturbation      [in percent]
8      Column wise values begin here.                   [] 
C 1    Axial z value from DIF3D                         [cm] 
C 2    Power for that AREA                              [:math:`Watt/cm^2`]
C 3    Fractional power for that AREA                   [normalized by total]
C 4    Doppler reactivity worth                         [reactivity/:math:`cm^2`]
C 5    Density reactivity worth                         [reactivity/:math:`cm^2`]   
====== ================================================ =========================
