MODULE cubes
  USE kinds, ONLY: wp => dp
  IMPLICIT NONE

  TYPE, PUBLIC :: cube
    PRIVATE
    CHARACTER (LEN=72) :: str1
    CHARACTER (LEN=72) :: str2
    REAL (KIND=wp) :: x_min, y_min, z_min, dx, dy, dz
    INTEGER :: n_x, n_y, n_z, N_atoms
    INTEGER, DIMENSION(:), POINTER :: z_atoms
    REAL (KIND=wp), DIMENSION(:), POINTER :: chrg, x, y, z
    REAL (KIND=wp), DIMENSION(:), POINTER :: array
    REAL (KIND=wp) :: dummy1, dummy2
  END TYPE cube

  PUBLIC :: cube_get, &
            cube_add, &
            cube_sub, &
            cube_int, &
            cube_cdz, &
            cube_del

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE cube_add
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE cube_sub
  END INTERFACE

  CONTAINS

  ! get input from target .cube file 
  SUBROUTINE cube_get (mycube, infile)
    CHARACTER(LEN=*), INTENT(IN) :: infile
    TYPE (cube), INTENT(OUT) :: mycube
    INTEGER :: i
    
    OPEN (UNIT=11, FILE=infile, STATUS="old", ACTION="read")           
                                                                      
    READ (UNIT=11, FMT=*) mycube%str1
    READ (UNIT=11, FMT=*) mycube%str2
    READ (UNIT=11, FMT=*) mycube%N_atoms, mycube%x_min, mycube%y_min, mycube%z_min
    READ (UNIT=11, FMT=*) mycube%n_x, mycube%dx, mycube%dummy1, mycube%dummy2
    READ (UNIT=11, FMT=*) mycube%n_y, mycube%dummy1, mycube%dy, mycube%dummy2
    READ (UNIT=11, FMT=*) mycube%n_z, mycube%dummy1, mycube%dummy2, mycube%dy

    ! allocate to each pointer variable the presumed size (amount of atoms)
    ALLOCATE (mycube%z_atoms(mycube%N_atoms), mycube%chrg(mycube%N_atoms))
    ALLOCATE (mycube%x(mycube%N_atoms), mycube%y(mycube%N_atoms), mycube%z(mycube%N_atoms))
    ALLOCATE (mycube%array(mycube%n_x*mycube%n_y*mycube%n_z))          

    ! loop over the lines containing Z, a charge, and x/y/z coords
    DO i=1, mycube%N_atoms
      READ(UNIT=11, FMT=*) mycube%z_atoms(i), mycube%chrg(i), mycube%x(i), mycube%y(i), mycube%z(i)
    END DO
                                     
    READ (UNIT=11, FMT=*) mycube%array

    CLOSE (11)                                                         
    
  END SUBROUTINE cube_get

  ! operation of adding together two cubes
  FUNCTION cube_add (mycube1, mycube2)
    TYPE(cube) :: cube_add
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
    

    cube_add%array = mycube1%array + mycube2%array   
  END FUNCTION cube_add

  ! operation of subtracting between two cubes
  FUNCTION cube_sub (mycube1, mycube2)
    TYPE(cube) :: cube_sub
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
    ! ...
  END FUNCTION cube_sub

  ! integration function
  FUNCTION cube_int (mycube)
    REAL (KIND=wp) :: cube_int
    TYPE (cube), INTENT(IN) :: mycube
    ! ...
  END FUNCTION cube_int

  ! "destroy" the cube and deallocate the memory
  SUBROUTINE cube_del (mycube)
    TYPE (cube), INTENT(IN) :: mycube
    ! ...
  END SUBROUTINE cube_del

  ! cdz
  SUBROUTINE cube_cdz (mycube, cdz)
    TYPE (cube), INTENT(IN) :: mycube
    REAL (KIND=wp) :: cdz
    ! ...
  END SUBROUTINE cube_cdz


END MODULE cubes