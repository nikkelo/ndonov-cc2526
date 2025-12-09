PROGRAM cda
  USE kinds, ONLY: wp => dp
  USE cubes
  IMPLICIT NONE

  TYPE (cube) :: rho_a, rho_b, rho_ab, rho_ref, drho
  REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: cdz

  ! using the cube_get subroutine, pull the data from the cube files
  CALL cube_get(rho_a, "../test/CuCO+/a.cube")
  CALL cube_get(rho_b, "../test/CuCO+/b.cube")
  CALL cube_get(rho_ab, "../test/CuCO+/ab.cube")

  ! calculate charge displacement according to formulas
  rho_ref = rho_a + rho_b 
  drho = rho_ab - rho_ref

  CALL cube_cdz(drho, cdz) 
  
  PRINT *, cdz

  ! delete the cubes from memory in case you want to start over
  CALL cube_del(rho_a)
  CALL cube_del(rho_b)
  CALL cube_del(rho_ab)

END PROGRAM cda