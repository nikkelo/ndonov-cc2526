PROGRAM verlet
  USE generic_energies
  USE kinds
  IMPLICIT NONE
  
  REAL (KIND = wp) :: tau, target_time, pos_diff
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: f_a, f_ap
  REAL (KIND = wp) :: m
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: pos_k, pos_kp, pos_temp
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: v_k, v_kp
  REAL (KIND = wp) :: E_k, E_p, E
  INTEGER :: k, n
  LOGICAL :: convergence

  ALLOCATE(f_a(3))
  ALLOCATE(f_ap(3))
  ALLOCATE(pos_k(3))
  ALLOCATE(pos_kp(3))
  ALLOCATE(pos_temp(3))
  ALLOCATE(v_k(3))
  ALLOCATE(v_kp(3))

  tau = 2_wp ! /s
  target_time = 120.0_wp ! Simulating for 2 minute
  f_a = [0.0_wp, 0.1_wp, 0.0_wp] ! /kg m s^-2 ; constant force with x, y and z components
  m = 1_wp ! /kg ; mass of particle
  pos_temp = [0.0_wp, 0.0_wp, 0.0_wp]
  n = 100 ! 10000 iteriations of steps k
  convergence = .FALSE. 
   
  DO WHILE (.NOT. convergence)
    pos_k = [0.0_wp, 0.0_wp, 0.0_wp]
    v_k = [0.0_wp, 0.0_wp, 0.0_wp] ! /m s^-1 ; speed of particle
    n = INT(target_time / tau)

    DO k = 1, n
      ! Compute positions, update force, compute velocities
      pos_kp = pos_k + tau*v_k + (tau**2.0_wp)*(f_a)/(2.0_wp*m)
      f_ap = f_a
      v_kp = v_k + tau/(2.0_wp*m)*(f_a+f_ap)

      PRINT *, "PARTICLE POSITIONS AT ITERATION", k, "AND TIMESTEP", tau*k, ":"
      PRINT *, pos_kp

      PRINT *

      ! Compute the kinetic, potential and total energy
      CALL kinetic_energy(m, pos_k(1), pos_k(2), pos_k(3), E_k)
      E_p = -DOT_PRODUCT(f_ap, pos_kp) ! potential energy calculation
      E = E_k + E_p

      PRINT *, "KINETIC, POTENTIAL AND TOTAL ENERGY OF THE SYSTEM:"
      PRINT *, E_k, E_p, E

      PRINT *
      PRINT *

      ! Compute Euclidean distance from the previous loop
      pos_diff = SQRT(SUM((pos_kp - pos_temp)**2))

      ! Update the positions and velocities
      pos_k = pos_kp
      v_k = v_kp

    ENDDO

    pos_temp = pos_k
    
    IF (pos_diff <= 0.01_wp .AND. tau < 0.2_wp) THEN ! Check if convergence has occured within a difference of 1cm
      convergence = .TRUE.
      PRINT *, "POSITIONS CONVERGED WITH 1 CM PRECISION AFTER 2 MINUTES, AT TIME STEP", tau
    ELSE ! If not converged, half the timestep tau and redo
      pos_temp = pos_kp
      tau = tau / 2.0_wp
    ENDIF
  ENDDO
 
  DEALLOCATE(f_a, f_ap, pos_k, pos_kp, v_k, v_kp, pos_temp)
END PROGRAM verlet
