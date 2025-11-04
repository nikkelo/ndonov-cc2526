PROGRAM verlet_harmonic
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)
  REAL (KIND = wp) :: tau
  REAL (KIND = wp) :: f_a, f_ap
  REAL (KIND = wp) :: m
  REAL (KIND = wp) :: x_k, x_kp, x_eq ! notation change from pos_k to x_k since we have only a single coordinate 
  REAL (KIND = wp) :: v_k, v_kp
  REAL (KIND = wp) :: k_spring
  INTEGER :: k, n ! keep in mind k is not the spring constant, just the step counter

  tau = 0.05_wp ! /s
  f_a = 0.1_wp ! /kg m s^-2 ; constant force with x, y and z components
  v_k = 0.0_wp ! /m s^-1 ; speed of particle
  m = 1_wp ! /kg ; mass of particle
  x_k = 1.0_wp ! position of x, displaced from equilibrium
  x_eq = 0.0_wp ! position of x at equilibrium
  n = 50 ! 50 iterations of steps k
  k_spring = 100.0_wp ! /kg s^2 ; spring constant 

  DO k = 1, n
    f_a = -k_spring * (x_k - x_eq)
    x_kp = x_k + tau*v_k + (tau**2)*(f_a)/(2.0_wp*m)
    f_ap = -k_spring * (x_kp - x_eq)
    v_kp = v_k + tau/(2.0_wp*m)*(f_a+f_ap)

    PRINT *, "PARTICLE POSITIONS AT ITERATION", k, "AND TIMESTEP", tau*k, ":"
    PRINT *, x_k

    PRINT *
    PRINT *

    x_k = x_kp
    v_k = v_kp

  ENDDO

  DEALLOCATE(f_a, f_ap, pos_k, pos_kp, v_k, v_kp)

END PROGRAM verlet_harmonic
