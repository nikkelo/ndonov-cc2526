PROGRAM verlet_LJ
  USE kinds, ONLY: wp => dp
  IMPLICIT NONE
  
  ! define parameters
  REAL (KIND = wp) :: tau  ! timestep, loop count and total number of steps
  INTEGER :: k, n ! loop parameters
  REAL (KIND = wp) :: r ! distance
  REAL (KIND = wp) :: sigma, epsilon, m ! parameters for solving the equation, read from the XYZ file 
  REAL (KIND = wp), DIMENSION(2,3) :: pos, pos_nxt, v, v_nxt, f, f_nxt ! position at current iteration, velocity at current iteration, force at current iteration, force at next iteration

  ! loop parameters
  tau = 1.0_wp
  k = 1
  n = 6000 

  ! given parameters by .xyz file
  sigma = 5.2186_wp
  epsilon = 0.000112991_wp

  ! assign mass, starting positions, velocities and forces by .xyz file. index(i, j) correpsonds to the i-th atom and its x/y/z component, respectively
  m = 20.1797_wp

  pos(1,:) = (/ 0.0_wp, 0.0_wp, 8.0_wp /) ! positions of the first Ne atom
  pos(2,:) = (/ 0.0_wp, 0.0_wp, 0.0_wp /) ! positions of the second Ne atom

  v(1,:) = (/ 0.0_wp, 0.0_wp, 0.0_wp /) ! velocities of the first Ne atom
  v(2,:) = (/ 0.0_wp, 0.0_wp, 0.0_wp /) ! velocities of the first Ne atom

  f(1,:) = (/ 0.0_wp, 0.0_wp, 0.0_wp /) ! forces of the first Ne atom
  f(2,:) = (/ 0.0_wp, 0.0_wp, 0.0_wp /) ! forces of the second Ne atom

  ! open neon.xyz for outputs
  OPEN (UNIT=12, FILE="neon.xyz", STATUS="replace", ACTION="write")

  DO k = 1, n

    ! calculate force at current positions
    r = SQRT(SUM((pos(1,:) - pos(2,:))**2.0_wp))
    f(1,:) = - ((((pos(1,:) - pos(2,:))/r)* 4.0_wp * epsilon * (-12.0_wp * (sigma**12)/(r**12) + 6.0_wp * (sigma**6)/(r**7))))
    f(2,:) = - ((((pos(2,:) - pos(1,:))/r)* 4.0_wp * epsilon * (-12.0_wp * (sigma**12)/(r**12) + 6.0_wp * (sigma**6)/(r**7))))

    ! update the positions of both atoms
    pos_nxt(1,:) = pos(1,:) + tau*v(1,:) + (tau**2)*f(1,:)/(2.0_wp*m)
    pos_nxt(2,:) = pos(2,:) + tau*v(2,:) + (tau**2)*f(2,:)/(2.0_wp*m)

    ! calculate force at new positions
    r = SQRT(SUM((pos(1,:) - pos(2,:))**2))
    f_nxt(1,:) = -((((pos_nxt(1,:) - pos_nxt(2,:))/r)* 4.0_wp * epsilon * (-12.0_wp * (sigma**12)/(r**12) + 6.0_wp * (sigma**6)/(r**7))))
    f_nxt(2,:) = -((((pos_nxt(2,:) - pos_nxt(1,:))/r)* 4.0_wp * epsilon * (-12.0_wp * (sigma**12)/(r**12) + 6.0_wp * (sigma**6)/(r**7))))

    ! update the velocities of both atoms
    v_nxt(1,:) = v(1,:) + tau/(2.0_wp*m)*(f(1,:) + f_nxt(1,:))
    v_nxt(2,:) = v(2,:) + tau/(2.0_wp*m)*(f(2,:) + f_nxt(2,:))

    ! print the current positions of the file 
    PRINT *, "PARTICLE POSITIONS AT ITERATION", k, "AND TIMESTEP", tau*k, ":"
      PRINT *, pos_nxt

    ! write to .xyz file output

    WRITE (UNIT=12, FMT=*) "2"
    WRITE (UNIT=12, FMT=*) "TIMESTEP: ", tau 
    WRITE (UNIT=12, FMT=*) "Ne", pos_nxt(1,:)
    WRITE (UNIT=12, FMT=*) "Ne", pos_nxt(2,:)

    ! update for next iteration
    pos = pos_nxt
    v = v_nxt

  ENDDO

  CLOSE (12)
      
END PROGRAM verlet_LJ
