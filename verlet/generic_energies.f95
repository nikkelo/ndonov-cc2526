MODULE generic_energies
   IMPLICIT NONE
   INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)

   PUBLIC :: kinetic_energy, potential_energy, total_energy

   CONTAINS 

   SUBROUTINE kinetic_energy(m, x, y, z, E_k)

      REAL(KIND=wp), INTENT(IN)  :: m, x, y, z   ! Inputs for kinetic energy formula
      REAL(KIND=wp), INTENT(OUT) :: E_k           ! Kinetic energy output

      E_k = 0.5_wp * m * (x**2 + y**2 + z**2)

   END SUBROUTINE kinetic_energy

   SUBROUTINE potential_energy(f, pos, E_p)

      REAL(KIND=wp), DIMENSION(3), INTENT(IN) :: f, pos ! Inputs for potential energy formula
      REAL(KIND=wp), INTENT(OUT) :: E_p ! Potential energy output
      E_p = -DOT_PRODUCT(f, pos)

   END SUBROUTINE potential_energy

   SUBROUTINE total_energy(K, P, E)

      REAL(KIND=wp), INTENT(IN) :: K, P ! Inputs for total energy calculation
      REAL(KIND=wp), INTENT(OUT) :: E ! Output for total energy 
      E = K + P

   END SUBROUTINE total_energy

END MODULE generic_energies