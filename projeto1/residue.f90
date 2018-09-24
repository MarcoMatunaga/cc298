subroutine residue(r_i,r_j)
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
if (which_diss == 1) then 
    !
    !
    D4_ksi(r_i,r_j,1) = dis_ksi4(r_i,r_j,1)
    D4_ksi(r_i,r_j,2) = dis_ksi4(r_i,r_j,2)
    D4_ksi(r_i,r_j,3) = dis_ksi4(r_i,r_j,3)
    D4_ksi(r_i,r_j,4) = dis_ksi4(r_i,r_j,4)
    !
    !
    D4_eta(r_i,r_j,1) = dis_eta4(r_i,r_j,1)
    D4_eta(r_i,r_j,2) = dis_eta4(r_i,r_j,2)
    D4_eta(r_i,r_j,3) = dis_eta4(r_i,r_j,3)
    D4_eta(r_i,r_j,4) = dis_eta4(r_i,r_j,4)
    !
    !
else if (which_diss == 2) then 
    !
    !
    D4_ksi(r_i,r_j,1) = dis_ksi2(r_i,r_j,1)
    D4_ksi(r_i,r_j,2) = dis_ksi2(r_i,r_j,2)
    D4_ksi(r_i,r_j,3) = dis_ksi2(r_i,r_j,3)
    D4_ksi(r_i,r_j,4) = dis_ksi2(r_i,r_j,4)
    !
    !
    D4_eta(r_i,r_j,1) = dis_eta2(r_i,r_j,1)
    D4_eta(r_i,r_j,2) = dis_eta2(r_i,r_j,2)
    D4_eta(r_i,r_j,3) = dis_eta2(r_i,r_j,3)
    D4_eta(r_i,r_j,4) = dis_eta2(r_i,r_j,4)
    !
    !
end if
!
!
if (r_i == 1 .or. r_j == 1) then 
                    !
                    !
residue1(r_i,r_j) = 0.50d0*( E_barra(r_i+1,r_j,1) - E_barra(r_i,r_j,1) ) &
                    + 0.50d0*( F_barra(r_i,r_j+1,1) - F_barra(r_i,r_j,1) )
                    !
                    !
residue2(r_i,r_j) = 0.50d0*( E_barra(r_i+1,r_j,2) - E_barra(r_i,r_j,2) ) &
                    + 0.50d0*( F_barra(r_i,r_j+1,2) - F_barra(r_i,r_j,2) )
                    !
                    !
residue3(r_i,r_j) = 0.50d0*( E_barra(r_i+1,r_j,3) - E_barra(r_i,r_j,3) ) &
                    + 0.50d0*( F_barra(r_i,r_j+1,3) - F_barra(r_i,r_j,3) )
                    !
                    !
residue4(r_i,r_j) = 0.50d0*( E_barra(r_i+1,r_j,4) - E_barra(r_i,r_j,4) ) &
                    + 0.50d0*( F_barra(r_i,r_j+1,4) - F_barra(r_i,r_j,4) )
                    !
                    !
else if (r_j == jmax .or. r_i == imax) then 
                    !
                    !
residue1(r_i,r_j) = 0.50d0*( E_barra(r_i,r_j,1) - E_barra(r_i-1,r_j,1) ) &
                    + 0.50d0*( F_barra(r_i,r_j,1) - F_barra(r_i,r_j-1,1) )
                    !
                    !
residue2(r_i,r_j) = 0.50d0*( E_barra(r_i,r_j,2) - E_barra(r_i-1,r_j,2) ) &
                    + 0.50d0*( F_barra(r_i,r_j,2) - F_barra(r_i,r_j-1,2) )
                    !
                    !
residue3(r_i,r_j) = 0.50d0*( E_barra(r_i,r_j,3) - E_barra(r_i-1,r_j,3) ) &
                    + 0.50d0*( F_barra(r_i,r_j,3) - F_barra(r_i,r_j-1,3) )
                    !
                    !
residue4(r_i,r_j) = 0.50d0*( E_barra(r_i,r_j,4) - E_barra(r_i-1,r_j,4) ) &
                    + 0.50d0*( F_barra(r_i,r_j,4) - F_barra(r_i,r_j-1,4) )
                    !
                    !   
end if
!
!
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
        if ( log10(abs(residue1(r_i,r_j))) > max_residue ) max_residue = log10(abs(residue1(r_i,r_j)))
        if ( log10(abs(residue2(r_i,r_j))) > max_residue ) max_residue = log10(abs(residue2(r_i,r_j)))
        if ( log10(abs(residue3(r_i,r_j))) > max_residue ) max_residue = log10(abs(residue3(r_i,r_j)))
        if ( log10(abs(residue4(r_i,r_j))) > max_residue ) max_residue = log10(abs(residue4(r_i,r_j)))
        !
        !
end subroutine residue