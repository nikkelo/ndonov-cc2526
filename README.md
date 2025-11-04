# ndonov-cc2526
A repository for the **Coding in Chemistry** course, as part of CDDC degre in the **University of Padova**, academic year 2025/2026.  

Contains primarily exercises from lessons as well as some additional experiments I do to better understand the material. I try to document everything I do within this README.md file, plus some extra comments within the code itself.

The course is taught by prof. Sergio Rampino, and there is a repository containing a lot of the instructions for the course: [coding-in-chemistry](https://github.com/srampinogroup/coding-in-chemistry).
# verlet

The files in the folder titled **verlet** contain Fortran files implementing the [Verlet algorithm](https://en.wikipedia.org/wiki/Verlet_integration), for the calculation of **particle trajectories/positions** in molecular dynamics symulations.

There are multiple files in the folder each with its own variation of the implementation:
- **verlet.f95**, which is the main file contains the simplest implementation with the force remaining constant throughout the simulation
- **verlet_harmonic.f95**, which models the moving particle as a harmonic oscillator
