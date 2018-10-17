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
max_residue = -100.0d0
!
if (which_diss == 1) then 
    !
    !
    D4_ksi(r_i,r_j,1) = dis_ksi4(r_i,r_j,1,eps_dis_e)
    D4_ksi(r_i,r_j,2) = dis_ksi4(r_i,r_j,2,eps_dis_e)
    D4_ksi(r_i,r_j,3) = dis_ksi4(r_i,r_j,3,eps_dis_e)
    D4_ksi(r_i,r_j,4) = dis_ksi4(r_i,r_j,4,eps_dis_e)
    !
    !
    D4_eta(r_i,r_j,1) = dis_eta4(r_i,r_j,1,eps_dis_e)
    D4_eta(r_i,r_j,2) = dis_eta4(r_i,r_j,2,eps_dis_e)
    D4_eta(r_i,r_j,3) = dis_eta4(r_i,r_j,3,eps_dis_e)
    D4_eta(r_i,r_j,4) = dis_eta4(r_i,r_j,4,eps_dis_e)
    !
    !
else if (which_diss == 2) then 
    !
    !
    D4_ksi(r_i,r_j,1) = dis_ksi2(r_i,r_j,1,eps_dis_e)
    D4_ksi(r_i,r_j,2) = dis_ksi2(r_i,r_j,2,eps_dis_e)
    D4_ksi(r_i,r_j,3) = dis_ksi2(r_i,r_j,3,eps_dis_e)
    D4_ksi(r_i,r_j,4) = dis_ksi2(r_i,r_j,4,eps_dis_e)
    !
    !
    D4_eta(r_i,r_j,1) = dis_eta2(r_i,r_j,1,eps_dis_e)
    D4_eta(r_i,r_j,2) = dis_eta2(r_i,r_j,2,eps_dis_e)
    D4_eta(r_i,r_j,3) = dis_eta2(r_i,r_j,3,eps_dis_e)
    D4_eta(r_i,r_j,4) = dis_eta2(r_i,r_j,4,eps_dis_e)
    !
    !
else if (which_diss == 3) then 
    !
    !
    call allocate_vars_non_linear
    call non_linear_dissipation
    !
    !
    D4_ksi(r_i,r_j,1) = non_lin_dis_ksi(r_i,r_j,eps2_ksi(r_i,r_j),eps4_ksi(r_i,r_j),1)
    D4_ksi(r_i,r_j,2) = non_lin_dis_ksi(r_i,r_j,eps2_ksi(r_i,r_j),eps4_ksi(r_i,r_j),2)
    D4_ksi(r_i,r_j,3) = non_lin_dis_ksi(r_i,r_j,eps2_ksi(r_i,r_j),eps4_ksi(r_i,r_j),3)
    D4_ksi(r_i,r_j,4) = non_lin_dis_ksi(r_i,r_j,eps2_ksi(r_i,r_j),eps4_ksi(r_i,r_j),4)
    !
    !
    D4_eta(r_i,r_j,1) = non_lin_dis_eta(r_i,r_j,eps2_eta(r_i,r_j),eps4_eta(r_i,r_j),1)
    D4_eta(r_i,r_j,2) = non_lin_dis_eta(r_i,r_j,eps2_eta(r_i,r_j),eps4_eta(r_i,r_j),2)
    D4_eta(r_i,r_j,3) = non_lin_dis_eta(r_i,r_j,eps2_eta(r_i,r_j),eps4_eta(r_i,r_j),3)
    D4_eta(r_i,r_j,4) = non_lin_dis_eta(r_i,r_j,eps2_eta(r_i,r_j),eps4_eta(r_i,r_j),4)
    !
    !
    call deallocate_vars_non_linear
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
        ! ***********************************
        ! posso juntar os dois blocos seguintes
        ! ***********************************
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
        if ( (abs(residue1(r_i,r_j))) > max_residue ) max_residue = (abs(residue1(r_i,r_j)))
        if ( (abs(residue2(r_i,r_j))) > max_residue ) max_residue = (abs(residue2(r_i,r_j)))
        if ( (abs(residue3(r_i,r_j))) > max_residue ) max_residue = (abs(residue3(r_i,r_j)))
        if ( (abs(residue4(r_i,r_j))) > max_residue ) max_residue = (abs(residue4(r_i,r_j)))
        !
        max_residue = log10(max_residue)
        !
        !
end subroutine residue