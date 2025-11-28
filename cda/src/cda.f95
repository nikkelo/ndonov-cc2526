PROGRAM cda
  USE kinds, ONLY: wp => dp
  USE cubes
  IMPLICIT NONE

  TYPE (cube) :: rho_a, rho_b, rho_ab, rho_ref, drho
  REAL (KIND=wp) :: cdz

  ! using the cube_get subroutine, pull the data from the cube files
  CALL cube_get(rho_a, "../test/CuCO+/a.cube")
  CALL cube_get(rho_b, "../test/CuCO+/b.cube")
  CALL cube_get(rho_ab, "../test/CuCO+/ab.cube")

  ! calculate charge displacement according to formulas
  


END PROGRAM cda
