# ndonov-cc2526
A repository for the **Coding in Chemistry** course, as part of CDDC degre in the **University of Padova**, academic year 2025/2026.  

Contains primarily exercises from lessons as well as some additional experiments I do to better understand the material. I try to document everything I do within this README.md file, plus some extra comments within the code itself.

The course is taught by prof. Sergio Rampino, and there is a repository containing a lot of the instructions for the course: [coding-in-chemistry](https://github.com/srampinogroup/coding-in-chemistry).
# verlet

The files in the folder titled **verlet** contain Fortran files implementing the [Verlet algorithm](https://en.wikipedia.org/wiki/Verlet_integration), for the calculation of **particle trajectories/positions** in molecular dynamics symulations.

There are multiple files in the folder each with its own variation of the implementation, or an additional subroutine to be used in the calculations for the individual variations. 

## Simplest Implementations
- **verlet.f95**, which is the main file contains the simplest implementation with the force remaining constant throughout the simulation. 
- **verlet_harmonic.f95**, which models the moving particle as a harmonic oscillator

Some additional lines were added to do exercises about checking the convergence at given points, calculating the energies at every point, doing some of the calculations from modules rather than in the main code, etc. 

## Lennard-Jones Potentials
- **verlet_LJ.f95**, which models atoms in 3D space interacting through pairwise additive Lennard-Jones potentials. The results are written into the "neon.xyz" file in XYZ format and can from there be used in VMD to visualize.  

## Supplementary Files
- **kinds.f95**, a module containing the two working precisions used (double precision and single precision).
- **generic_energies.f95**, another module containing the energy calculations for the verlet.f95 implementation
- **neon.xyz**, contains the output from the interacting Neon atoms for the LJ implementation. Parameters sigma and epsilon sourced from [this paper](doi.org/10.1063/1.4796144)

# cda

The files in the folder titled **cda** contain Fortran files for implementing a Charge-Displacement Analysis on a set of test .CUBE files referring to CuCO+. 

The main file located at **src/cda.f95** contains a compact implementation, utilizing the subroutines found in the **src/cubes.f95** module.

The following subroutes can be found in **src/cubes.f95** which are later used to compute the charge displacement in the main file:
- cube_get: read a .cube file
- cube_add: operation of addition for cubes
- cube_sub: operation of subtraction for cubes
- cube_unwrap: "unwrap" a cube into a 3D array from the 1D charge array found at the end of the file
- cube_int: integrate over xy
- cube_cdz: compute the charge displacement along z
- cube_del: unallocate a cube from memory

# py

The files in the folder titled **py** are jupyter notebooks related to various procedures for the use of machine learning in chemistry (scikit-learn).

