PROGRAM verlet
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)
  REAL (KIND = wp) :: tau
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: f_a, f_ap
  REAL (KIND = wp) :: m
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: pos_k, pos_kp
  REAL (KIND = wp), DIMENSION(:), ALLOCATABLE :: v_k, v_kp
  REAL (KIND = wp) :: E_k, E_p, E
  INTEGER :: k, n

  ALLOCATE(f_a(3))
  ALLOCATE(f_ap(3))
  ALLOCATE(pos_k(3))
  ALLOCATE(pos_kp(3))
  ALLOCATE(v_k(3))
  ALLOCATE(v_kp(3))

  tau = 0.02_wp ! /s
  f_a = [0.0_wp, 0.1_wp, 0.0_wp] ! /kg m s^-2 ; constant force with x, y and z components
  v_k = [0.0_wp, 0.0_wp, 0.0_wp] ! /m s^-1 ; speed of particle
  m = 1_wp ! /kg ; mass of particle
  pos_k = [0.0_wp, 0.0_wp, 0.0_wp]
  n = 1000 ! 10000 iterations of steps k

  DO k = 1, n
    pos_kp = pos_k + tau*v_k + (tau**2)*(f_a)/(2.0_wp*m)
    f_ap = f_a
    v_kp = v_k + tau/(2.0_wp*m)*(f_a+f_ap)

    E_k = 0.5_wp*m*(v_kp(1)**2 + v_kp(2)**2 + v_kp(3)**2) ! kinetic energy calculation
    E_p = -DOT_PRODUCT(f_ap, pos_kp) ! potential energy calculation
    E = E_k + E_p


    PRINT *, "PARTICLE POSITIONS AT ITERATION", k, "AND TIMESTEP", tau*k, ":"
    PRINT *, pos_kp

    PRINT *

    PRINT *, "KINETIC, POTENTIAL AND TOTAL ENERGY OF THE SYSTEM:"
    PRINT *, E_k, E_p, E

    PRINT *
    PRINT *

    pos_k = pos_kp
    v_k = v_kp

  ENDDO

  DEALLOCATE(f_a, f_ap, pos_k, pos_kp, v_k, v_kp)
END PROGRAM verlet
