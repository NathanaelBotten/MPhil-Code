# Dirac Solver

## Solving the Dirac Equation in a Carbon Nucleus

### Nathanael Botten, 2023; Modified from the work of Dr. Jon Carroll 2010, revised in 2023.

This code employs the Runge-Kutta method in order to solve the Dirac equation for the nuclear environment, as described in the QMC model.

## Running the program

1. Download or clone this repository
1. Edit the `Makefile` to use whichever Fortran compiler you have. This has been tested with `gfortran`
1. Ensure a clean start with `make clean` which removes all `.o` and `.mod` files as well as the compiled binary
1. Compile the code with `make` (which automatically builds `make Dirac`)
1. Run the program with `./Dirac` which will perform the integration and produce the `Wavefunctions.dat` data file. 
   This will take some time, but provides progress updates to the terminal
1. Use your favourite plotting program to visualise the `F` and `G` components of the spinor
1. The code will return the value of the 3d-2p states to the terminal

## Editting the Code

Included in the file FindEigenvalue.f90 is a parameter in the subroutine `myDerivs` that switches between investigating the point Coulomb, finite Coulomb, and nuclear potentials. 
The form for these potentials are included in the file funcs.f90.

### References
Carroll, Jonathan David; Thomas, Anthony William; Rafelski, J.; Miller, G. A.
Nonperturbative relativistic calculation of the muonic hydrogen spectrum Physical Review
A, 2011; 84(1):012506 
