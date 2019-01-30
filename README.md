# MSR1D
Modified point kinetics solver for MSR reactor transient calculations.  

## Building
To build from the source code all that is needed is a fortran compiler. For instance [gfortran](https://gcc.gnu.org/wiki/GFortran)
A makefile is used to build the executable.  The makefile assumes you have gfortran in your path
```
$ cd [src](src)
$ make
```
Note: Regression testing and visualization utilities requires python 2.7+ with matplotlib & numpy

## Documentation
ocumentation can be found: [Click here for documentation](https://msr-1d.readthedocs.io/en/latest)

## Testing
Basic regression tests can be found by going to:
```
$ cd [regression](regression) 
```
And run using the python script:
```
$ python BasicRegressionTests.py
```
