subroutine residue
    use vars
    implicit none
! *************************************
! it is possible to create a subroutine residue to calculate the RHS
! *************************************
!
integer r_i, r_j
! residue - pick the maximum residue of them 
! which means the norma_infinity
!
if (which_diss == 1) call artificial_dissipation_D4
if (which_diss == 2) call artificial_dissipation_D2

max_residue = -1.0d0
do r_j = 2, jmax - 1
        do r_i = 2, imax - 1
residue1(r_i,r_j) = 0.50d0*(E_barra(r_i+1,r_j,1) - E_barra(r_i-1,r_j,1))  & 
              + 0.50d0*(F_barra(r_i,r_j+1,1) - F_barra(r_i,r_j-1,1)) 
        !
        !
residue2(r_i,r_j) = 0.50d0*(E_barra(r_i+1,r_j,2) - E_barra(r_i-1,r_j,2))  &
              + 0.50d0*(F_barra(r_i,r_j+1,2) - F_barra(r_i,r_j-1,2)) 
        !
        !
residue3(r_i,r_j) = 0.50d0*(E_barra(r_i+1,r_j,3) - E_barra(r_i-1,r_j,3)) & 
              + 0.50d0*(F_barra(r_i,r_j+1,3) - F_barra(r_i,r_j-1,3)) 
        !
        !
residue4(r_i,r_j) = 0.50d0*(E_barra(r_i+1,r_j,4) - E_barra(r_i-1,r_j,4)) &
              + 0.50d0*(F_barra(r_i,r_j+1,4) - F_barra(r_i,r_j-1,4)) 
        !
        !
        residue1(r_i,r_j) = residue1(r_i,r_j) + metric_jacobian(r_i,r_j)*(D4_ksi(r_i,r_j,1) + D4_eta(r_i,r_j,1))
        residue2(r_i,r_j) = residue2(r_i,r_j) + metric_jacobian(r_i,r_j)*(D4_ksi(r_i,r_j,2) + D4_eta(r_i,r_j,2))
        residue3(r_i,r_j) = residue3(r_i,r_j) + metric_jacobian(r_i,r_j)*(D4_ksi(r_i,r_j,3) + D4_eta(r_i,r_j,3))
        residue4(r_i,r_j) = residue4(r_i,r_j) + metric_jacobian(r_i,r_j)*(D4_ksi(r_i,r_j,4) + D4_eta(r_i,r_j,4))
        !
        !
        residue1(r_i,r_j) = delta_t(r_i,r_j)*residue1(r_i,r_j)
        residue2(r_i,r_j) = delta_t(r_i,r_j)*residue2(r_i,r_j)
        residue3(r_i,r_j) = delta_t(r_i,r_j)*residue3(r_i,r_j)
        residue4(r_i,r_j) = delta_t(r_i,r_j)*residue4(r_i,r_j)
        !
        !
        !
        write(*,*) E_barra(r_i,r_j,1), F_barra(r_i,r_j,1)
        if ( abs(residue1(r_i,r_j)) > max_residue ) max_residue = log10(abs(residue1(r_i,r_j)))
        if ( abs(residue2(r_i,r_j)) > max_residue ) max_residue = log10(abs(residue2(r_i,r_j)))
        if ( abs(residue3(r_i,r_j)) > max_residue ) max_residue = log10(abs(residue3(r_i,r_j)))
        if ( abs(residue4(r_i,r_j)) > max_residue ) max_residue = log10(abs(residue4(r_i,r_j)))
        !
        !
        end do
end do
!
!
!
end subroutine residue