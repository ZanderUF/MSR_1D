.. _Input: 

Input File Guide
================

The structure of the input files is broken up into blocks with character sequences seperating the blocks
The first line in the file is a comment, this should be something meaningful to distinguish the problem from others

.. Parameter block::

Parameter Block
---------------

Begins with **read parm** and ends with **end parm**

=========  ===========  ================================================
Parameter  Value        Description
=========  ===========  ================================================
time       yes or no    Perform a time dependent calculation
step       yes or no    Perform a step perturbation or not
ramp       yes or no    Perform a ramp perturbation or not
zag        yes or no    Perform a zig-zag perturbation or not
del        Double       Indicate fixed time step size, ex: 1E-5
tmax       Real         Indicate the total length of the simulation
tin        Double       Beginning point of the perturbation 
nem        Integer      Number of equal spaced elements in the model
npe        Integer      Number of nodes per element 
pipe       Integer      Number of elements in the outer circuit
area       Real         Area of the fuel core [cm^2] 
mflow      Real         Mass flow rate [gm/cm^3]
tpow       Real         Total axial power [n/cm]
nitr       Real         Set number of maximum nonlinear iterations 
elem       Real         Element length [cm]
ndg        Integer      Number of delayed neutron groups
nmat       Integer      Number of fissional materials
gen        Real         Neutron generation time [s^-2]
reac       Real         Reactivity insertion for step and ramp
=========  ===========  ================================================

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
       mflow=000
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

