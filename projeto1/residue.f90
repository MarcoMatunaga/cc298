subroutine compute_residue(r_i,r_j)
    use vars
    use functions
    implicit none
!
!
integer r_i, r_j
!
! residue - pick the maximum residue of them 
! which means the norma_infinity
!
call artificial_dissipation(r_i,r_j)
max_residue = -100.0d0
residue = 0.0d0
!
!
residue(r_i,r_j,1) = 0.50d0*(E_barra(r_i+1,r_j,1) - E_barra(r_i-1,r_j,1)) & 
                    + 0.50d0*(F_barra(r_i,r_j+1,1) - F_barra(r_i,r_j-1,1)) 
        !
        !
residue(r_i,r_j,2) = 0.50d0*(E_barra(r_i+1,r_j,2) - E_barra(r_i-1,r_j,2)) &
                    + 0.50d0*(F_barra(r_i,r_j+1,2) - F_barra(r_i,r_j-1,2)) 
        !
        !
residue(r_i,r_j,3) = 0.50d0*(E_barra(r_i+1,r_j,3) - E_barra(r_i-1,r_j,3)) & 
                    + 0.50d0*(F_barra(r_i,r_j+1,3) - F_barra(r_i,r_j-1,3)) 
        !
        !
residue(r_i,r_j,4) = 0.50d0*(E_barra(r_i+1,r_j,4) - E_barra(r_i-1,r_j,4)) &
                    + 0.50d0*(F_barra(r_i,r_j+1,4) - F_barra(r_i,r_j-1,4)) 
        !
        !
        residue(r_i,r_j,1) = delta_t(r_i,r_j)*residue(r_i,r_j,1) + (D4_ksi(r_i,r_j,1) + D4_eta(r_i,r_j,1))
        residue(r_i,r_j,2) = delta_t(r_i,r_j)*residue(r_i,r_j,2) + (D4_ksi(r_i,r_j,2) + D4_eta(r_i,r_j,2))
        residue(r_i,r_j,3) = delta_t(r_i,r_j)*residue(r_i,r_j,3) + (D4_ksi(r_i,r_j,3) + D4_eta(r_i,r_j,3))
        residue(r_i,r_j,4) = delta_t(r_i,r_j)*residue(r_i,r_j,4) + (D4_ksi(r_i,r_j,4) + D4_eta(r_i,r_j,4))
        !
        !
        if ( (abs(residue(r_i,r_j,1))) > max_residue ) max_residue = (abs(residue(r_i,r_j,1)))
        if ( (abs(residue(r_i,r_j,2))) > max_residue ) max_residue = (abs(residue(r_i,r_j,2)))
        if ( (abs(residue(r_i,r_j,3))) > max_residue ) max_residue = (abs(residue(r_i,r_j,3)))
        if ( (abs(residue(r_i,r_j,4))) > max_residue ) max_residue = (abs(residue(r_i,r_j,4)))
        !
        max_residue = log10(max_residue)
        if (isnan(max_residue)) stop 'Divergiu :('
        !
        !
end subroutine compute_residue