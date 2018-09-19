.. _Installation:

Installation Guide
==================
Installation has only been tested on a Mac to date

--------------
Access to Code
--------------

All relevant code is located on Github and managed on there.  Simple go to

https://github.com/ZanderUF/MSR_1D

You can clone a version of the repo with::

    $ git clone https://github.com/ZanderUF/MSR_1D.git

------------------------
Install on Mac
------------------------

Installation has only been tested with GNU Fortran version 5.4.0
Requires LAPACK, tested with version 3.8.0
http://www.netlib.org/lapack/#_lapack_version_3_8_0_2

If you change the install location of LAPACK be sure to update the makefile in the MSR1D install directory.  Right now LAPACK is assumed to be two directories up from the *src* directory::
    
    ../../lapack-3.8.0

All of the source code is locate in::
    
    /MSR_1D/PrecursorSolve/src

A simple makefile is used to build an exectuable in the *src* directory. 

To install just use make in the *src* directory::
   
    $ make

Note, only tested with GNU Make 3.81

This will make an executable **a.out** by default

To clean up the installation directory simple execute the folling in the *src* directory::
    
    $ make clean

