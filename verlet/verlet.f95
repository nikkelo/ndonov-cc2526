PROGRAM verlet
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)
  REAL (KIND = wp) :: tau
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: f_a, f_ap
  REAL (KIND = wp) :: m
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: pos_k, pos_kp
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: v_k, v_kp
  REAL (KIND = wp) :: E_k
  INTEGER :: k, n

  ALLOCATE(f_a(3))
  ALLOCATE(pos_k(3))
  ALLOCATE(pos_kp(3))
  ALLOCATE(v_k(3))
  ALLOCATE(v_kp(3))

  tau = 0.002_wp ! /s
  f_a = [0.0_wp, 0.1_wp, 0.0_wp] ! /kg m s^-2 ; constant force with x, y and z components
  m = 1_wp ! /kg ; mass of particle
  pos_k = [0.0_wp, 0.0_wp, 0.0_wp]
  n = 100000 ! 100 iterations of steps k

  DO k = 1, n
    pos_kp = pos_k + tau*v_k + (tau**2)*(f_a)/(2.0_wp*m)
    f_ap = f_a
    v_kp = v_k + tau/(2.0_wp*m)*(f_a+f_ap)

    E_k = 0.5_wp*m*(v_k(1)**2 + v_k(2)**2 + v_k(3)**2)
    
    PRINT *, "The x, y and z positions of the particle at the", k, "-th iteration at time =", tau*k, "are:"
    PRINT *, pos_k

    PRINT *, "The potential, kinetic and total energies of the particle at this step are, respectively:"
    PRINT *, E_k

    pos_k = pos_kp
    v_k = v_kp

  ENDDO

END PROGRAM verlet
