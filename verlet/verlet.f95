PROGRAM verlet
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)
  REAL (KIND = wp) :: tau
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: f_a, f_ap
  REAL (KIND = wp) :: m
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: pos_k, pos_kp
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: v_k, v_kp
  INTEGER :: k, n

  ALLOCATE(f_a(3))
  ALLOCATE(pos_k(3))
  ALLOCATE(pos_kp(3))
  ALLOCATE(v_k(3))
  ALLOCATE(v_kp(3))

  tau = 0.2_wp ! /s
  f_a = [0.0, 0.1, 0.0] ! /kg m s^-2 ; constant force with x, y and z components
  m = 1_wp ! /kg ; mass of particle
  pos_k = [0.0, 0.0, 0.0]
  n = 10 ! 10 iterations of steps k

  DO k = 1, n
    pos_kp = pos_k + tau*v_k + (tau**2)*(f_a)/(2.0_wp*m)
    f_ap = f_a
    v_kp = v_k + tau/(2.0_wp*m)*(f_a+f_ap)

    PRINT *, "The x, y and z positions of the particle at the", k, "-th iteration at time =", tau*k, "are:"
    PRINT *, pos_k

    pos_k = pos_kp
    v_k = v_kp

  ENDDO

END PROGRAM verlet
