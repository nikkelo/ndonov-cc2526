PROGRAM verlet
  USE generic_energies
  USE kinds
  IMPLICIT NONE
  
  ! define parameters
  REAL (KIND = wp) :: tau ! /s, timestep
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: f, f_nxt ! /kg m s^-2, force vector (x, y, z) 
  REAL (KIND = wp) :: m ! /kg, mass of particle
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: pos, pos_nxt ! /m, position vector (x, y, z)
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: v, v_nxt ! /m s^-1, velocity vector (x, y, z)
  REAL (KIND = wp) :: E_k, E_p, E ! /J, kinetic, potential, and total energy
  INTEGER :: k, n ! loop parameters

  ! allocation of variables
  ALLOCATE(f(3))
  ALLOCATE(f_nxt(3))
  ALLOCATE(pos(3))
  ALLOCATE(pos_nxt(3))
  ALLOCATE(v(3))
  ALLOCATE(v_nxt(3))

  ! loop parameters
  tau = 2.0_wp
  n = 10 

  ! initial conditions 
  m = 1.0_wp 
  f = [0.0_wp, 0.1_wp, 0.0_wp]
  pos = [0.0_wp, 0.0_wp, 0.0_wp]
  v = [0.0_wp, 0.0_wp, 0.0_wp]

  ! main loop for calculating further iterations
  DO k = 1, n

    ! update positions
    pos_nxt = pos + tau*v + (tau**2.0_wp)*(f)/(2.0_wp*m)

    ! calculate force at new position
    f_nxt = f

    ! calculate velocity at new position
    v_nxt = v + tau/(2.0_wp*m)*(f + f_nxt)

    ! print current positions
    PRINT *, "PARTICLE POSITIONS AT ITERATION", k, "AND TIMESTEP", tau*k, ":"
    PRINT *, pos_nxt

    PRINT *

    ! compute the kinetic, potential and total energy
    CALL kinetic_energy(m, v(1), v(2), v(3), E_k)
    CALL potential_energy(f, pos, pos_nxt, E_p)
    CALL total_energy(E_k, E_p, E)

    PRINT *, "KINETIC, POTENTIAL AND TOTAL ENERGY OF THE SYSTEM:"
    PRINT *, E_k, E_p, E

    PRINT *
    PRINT *

    ! update position and velocity for next iteration
    pos = pos_nxt
    v = v_nxt

  ENDDO
 
  DEALLOCATE(f, f_nxt, pos, pos_nxt, v, v_nxt)

END PROGRAM verlet
