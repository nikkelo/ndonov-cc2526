# ndonov-cc2526
A repository for the **Coding in Chemistry** course, as part of CDDC degre in the **University of Padova**, academic year 2025/2026.  

Contains primarily exercises from lessons as well as some additional experiments I do to better understand the material. I try to document everything I do within this README.md file, plus some extra comments within the code itself.

The course is taught by prof. Sergio Rampino, and there is a repository containing a lot of the instructions for the course: [coding-in-chemistry](https://github.com/srampinogroup/coding-in-chemistry).
# verlet

The files in the folder titled **verlet** contain Fortran files implementing the [Verlet algorithm](https://en.wikipedia.org/wiki/Verlet_integration), for the calculation of **particle trajectories/positions** in molecular dynamics symulations.

The snippets below contain the code found in the loop that calculates the particle position at different iterations. 

## verlet.f95
```
DO k = 1, n

    pos_kp = pos_k + tau*v_k + (tau**2)*(f_a)/(2.0_wp*m)
    f_ap = f_a
    v_kp = v_k + tau/(2.0_wp*m)*(f_a+f_ap)

    pos_k = pos_kp
    v_k = v_kp

ENDDO
```
This is the simplest version of the algorithm assuming that $f^{(a,x)}_k = f^{(a,x)}_{k+1}$. 
