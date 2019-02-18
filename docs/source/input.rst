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
DebuggIt   yes or no    Turn on DEBUG printing - prints a lot of stuff
ReadDif3   yes or no    Read DIF3D external file for power profile and such 
NlIters    Integer      Set number of maximum nonlinear iterations 
Material   Integer      Number of fissional materials
NumDelay   Integer      Number of delayed neutron groups
TotalPow   Real         Total axial power [:math:`n/s-cm`]
=========  ===========  ==========================================================

.. Time block::

Time Block
----------

Begins with **read time** and ends with **end time**

=========  ===========  ==========================================================
Parameter  Value        Description
=========  ===========  ==========================================================
TDMethod   integer      0 - forward Euler 1 - backward Euler       
TimeSolv   yes or no    Perform a time dependent calculation
TimeStep   Double       Indicate fixed time step size, ex: 1E-5
EndTime    Real         Indicate the total length of the simulation
StrtTime   Double       Beginning point of the calculation 
SaveTime   Real         Interval to write out spatial solution files to
=========  ===========  ==========================================================

.. Perturbation block::

Perturbation Block
------------------

Begins with **read perturbation** and ends with **end perturbation**

=========  ===========  ==========================================================
Parameter  Value        Description
=========  ===========  ==========================================================
Feedback   Integer      0 - (default) no feedback 1 - temperature 2 - instant beta
StepPert   yes or no    Perform a step perturbation or not
StrtStep   Double       Start time for step perturbation
EndStep    Double       End time for step perturbation
RampPer    yes or no    Perform a ramp perturbation or not
StrtRamp   Double       Start time for ramp perturbation
EndRamp    Double       End time for ramp perturbation
ZaggPert   yes or no    Perform a zig-zag perturbation or not
Reactiv    Real         Reactivity insertion for step and ramp
TimeCons   Real         Time constant for mass flow rate reduction
PerFlow    Real         Percent reduction in flow rate over 1/TimeCons  
MassFlow   Real         Mass flow rate [:math:`g/cm^3`]
GenTime    Real         Neutron generation time [:math:`s^-2`]
=========  ===========  ==========================================================

.. Mesh block::

Mesh Block
------------------

Begins with **read mesh** and ends with **end mesh**

=========  ===========  ==========================================================
Parameter  Value        Description
=========  ===========  ==========================================================
ElemSize   Integer      Size of elements
NumElems   Integer      Number of equal spaced elements in the model
NumNodes   Integer      Number of nodes per element 
FuelInlt   Integer      Starting element of the inlet plenum
CoreStrt   Integer      Starting element of the main core region
CoreEnd    Integer      Final element of the main core region
FuelOutl   Integer      Final element of the outlet plenum
StartHex   Integer      Heat exchanger starting element location
EndHexch   Integer      Heat exchanger ending element location
CoreArea   Real         Area of the fuel core [:math:`cm^2`] 
PipArea    Real         Area of the piping [:math:`cm^2`]
HexcArea   Real         Area of the heat exchanger  [:math:`cm^2`]
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

    Heat Exchanger Overcooling Model of Chrloride type problem
    read parm
       DebuggIt=no
       ReadDif3=no
       NumDelay=6
       Material=1
       NlIters=50
       TotalPow=10
    end parm
    
    read time
       TimeSolv=yes
       TDMethod=1
       TimeStep=1E-3
       EndTime=20.0
       StrtTime=0.0 
       SaveTime=10.0
    end time
    
    read pert
       Feedback=0
       StepPert=yes
       RampPert=no
       ZaggPert=no
       StrtStep=0.0
       EndStep=1.0
       StrtRamp=0.0
       EndRamp=0.0
       Reactiv=0.003
       TimeCons=0.0
       PerFlow=0.0
       MassFlow=0.0
       GenTime=2.0E-5
    end pert
    
    read mesh
       ElemSize=1.0
       NumElems=10
       NumNodes=3
       FuelInlt=1
       CoreStrt=2
       CoreEnd=9
       FuelOutl=10
       StartHex=1
       EndHExch=1
       CoreArea=7.49E4
       PipeArea=7.49E4
       HexcArea=100000
    end mesh
    read delay
      mat=1
      alam=0.0127 0.0317 0.115 0.311 1.4 3.87 end
      beta=2.66E-4 1.491E-3 1.316E-3 2.849E-3 8.96E-4 1.82E-4  end
    end delay
    
.. Input File DIF3D Values::


Input file structure to read in from DIF3D
------------------------------------------

Values from DIF3D can be read in and projected onto the domain in the 1D problem
Right now it reads in the power, power fraction, doppler reactivity worth, and density reactivity worth. These are all spatially dependent.  These values are integrated across the AREA specificed in the DIF3D problem. 
It will read in the values and project them up to that axial height.  For example:

 =======  =======  =======   ========  ==========    =========   ==========
 Volume   Area     Z-Coord   Power     Frac Power    Doppler      Expansion
 =======  =======  =======   ========  ==========    =========   ==========
 200      20       10        0.00E+00  0.0+00        0            0   
 200      20       20        1.85E+06  6.1-03        -3.31E-07   -1.02E-06
 200      20       30        3.08E+06  1.0-02        -1.01E-06   -5.11E-07
 =======  =======  =======   ========  ==========    =========   ==========

This would project from locations 10 - 20 a power of 1.85E+06/(20-10).
Similalrly, from 20 - 30 a power of 3.08E+06/(30-20)

The file read in is assumed to be titled: dif3d_values.txt
It assumes the file is read in as follows:

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
8      Column wise values begin here.                    
C 1    Axial z value from DIF3D                         [cm] 
C 2    Volume                                           [:math:`cm^3`] 
C 3    Cross sectional area                             [:math:`cm^2`]
C 4    Power for that AREA                              [:math:`Watt/cm^2`]
C 5    Fractional power for that AREA                   [normalized by total]
C 6    Doppler reactivity worth                         [reactivity/:math:`cm^2`]
C 7    Density reactivity worth                         [reactivity/:math:`cm^2`]   
====== ================================================ =========================
