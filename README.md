# MSR1D
Modified point kinetics solver for MSR reactor transient calculations. Uses a 1D discontinous Galerkin approach for discretizing the temperature and fluid dyanmics equations.  

## Building
Requires LAPACK to be built.  See [lapack-3.8.0](lapack-3.8.0) for details on how to build
To build the source code all that is needed is a fortran compiler. For instance [gfortran](https://gcc.gnu.org/wiki/GFortran)
A makefile is used to build the executable.  The makefile assumes you have gfortran in your path
```
cd src
make
```
Note: Regression testing and visualization utilities requires python 2.7+ with matplotlib & numpy

## Documentation
Documentation maintained with Readthedocs. [Click here for documentation](https://msr-1d.readthedocs.io/en/latest)

## Testing
Basic regression tests can be found by going to:
```
cd regression
```
And run using the python script:
```
python BasicRegressionTests.py
```

## Utilities
The utilities folder has utilized for automatting steady state calculations to tabuluate various parameters.
Plotting of solution vectors is made possible with a series of pyhton scripts
