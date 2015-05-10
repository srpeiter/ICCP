# ICCP

ICCP is an International Course in Computational Physics
for physics Students from Delft University of technology. 

The first project is Molecular Dynamics. Here we simulated the behaviour of a noble gas.
The Codes for this project can be found in the folder Molecular dynamics in file ArgonSandbox_temp.
The Callaborator: Rick Vink

The second project is a quantum mechanics. Here we analyze the behaviour of 
bringing two hydrogen atoms close to each other. We solve the dissocaition energy for the hydrogen molecule with the quantum variatioanl method. This is done with a a paticular method of the Monte Carlo simulation also known as the metropolis walker. The algoritms are written in C++ and are located in the quantum mechanics folder (includes and sources).
Collaborator: Mario Gely

The third project is also a quantum mechanical simulation. Here we simulate the evolution of a quantum mechanical wavefunction. We solved the 1D and 2D time dependent Schrodinger equation numerically. For the 1D problem we used a Harmonic potential placed in the center of the domain. For the 2D case we did not use a potential barrier. Thus in the 2D case we only simulated the evolution of the wavefunction with the potential set to zero everywhere.  The algoritms are written in C++ and are located in the Quantum Dynamics folder. Type "make" in the working directory to run the simulation. Make sure your device support C++11 and gnuplot! Note that the libraries used here are: LAPACK and SUPERLU.
Collaborator: Mario Gely

