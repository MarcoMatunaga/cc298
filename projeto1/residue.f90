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
    if (r_i == 2 .or. r_i == imax - 1) then
            D4_ksi(r_i,r_j,1) = dis_ksi2(r_i,r_j,1,eps_dis_e)
            D4_ksi(r_i,r_j,2) = dis_ksi2(r_i,r_j,2,eps_dis_e)
            D4_ksi(r_i,r_j,3) = dis_ksi2(r_i,r_j,3,eps_dis_e)
            D4_ksi(r_i,r_j,4) = dis_ksi2(r_i,r_j,4,eps_dis_e)
    else if (r_j == 2 .or. r_j == jmax - 1) then
            D4_eta(r_i,r_j,1) = dis_eta2(r_i,r_j,1,eps_dis_e)
            D4_eta(r_i,r_j,2) = dis_eta2(r_i,r_j,2,eps_dis_e)
            D4_eta(r_i,r_j,3) = dis_eta2(r_i,r_j,3,eps_dis_e)
            D4_eta(r_i,r_j,4) = dis_eta2(r_i,r_j,4,eps_dis_e)
    else 
        Dnonlinear_plus   = non_lin_dis_ksi(r_i,r_j,eps2_ksi(r_i,r_j),eps4_ksi(r_i,r_j),1)
        Dnonlinear_minus  = non_lin_dis_ksi(r_i-1,r_j,eps2_ksi(r_i-1,r_j),eps4_ksi(r_i-1,r_j),1)
        D4_ksi(r_i,r_j,1) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_ksi(r_i,r_j,eps2_ksi(r_i,r_j),eps4_ksi(r_i,r_j),2)
        Dnonlinear_minus  = non_lin_dis_ksi(r_i-1,r_j,eps2_ksi(r_i-1,r_j),eps4_ksi(r_i-1,r_j),2)
        D4_ksi(r_i,r_j,2) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_ksi(r_i,r_j,eps2_ksi(r_i,r_j),eps4_ksi(r_i,r_j),3)
        Dnonlinear_minus  = non_lin_dis_ksi(r_i-1,r_j,eps2_ksi(r_i-1,r_j),eps4_ksi(r_i-1,r_j),3)
        D4_ksi(r_i,r_j,3) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_ksi(r_i,r_j,eps2_ksi(r_i,r_j),eps4_ksi(r_i,r_j),4)
        Dnonlinear_minus  = non_lin_dis_ksi(r_i-1,r_j,eps2_ksi(r_i-1,r_j),eps4_ksi(r_i-1,r_j),4)
        D4_ksi(r_i,r_j,4) = -Dnonlinear_plus + Dnonlinear_minus
        !
        !
        Dnonlinear_plus   = non_lin_dis_eta(r_i,r_j,eps2_eta(r_i,r_j),eps4_eta(r_i,r_j),1)
        Dnonlinear_minus  = non_lin_dis_eta(r_i,r_j-1,eps2_eta(r_i,r_j-1),eps4_eta(r_i,r_j-1),1)
        D4_eta(r_i,r_j,1) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_eta(r_i,r_j,eps2_eta(r_i,r_j),eps4_eta(r_i,r_j),2)
        Dnonlinear_minus  = non_lin_dis_eta(r_i,r_j-1,eps2_eta(r_i,r_j-1),eps4_eta(r_i,r_j-1),2)
        D4_eta(r_i,r_j,2) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_eta(r_i,r_j,eps2_eta(r_i,r_j),eps4_eta(r_i,r_j),3)
        Dnonlinear_minus  = non_lin_dis_eta(r_i,r_j-1,eps2_eta(r_i,r_j-1),eps4_eta(r_i,r_j-1),3)
        D4_eta(r_i,r_j,3) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_eta(r_i,r_j,eps2_eta(r_i,r_j),eps4_eta(r_i,r_j),4)
        Dnonlinear_minus  = non_lin_dis_eta(r_i,r_j-1,eps2_eta(r_i,r_j-1),eps4_eta(r_i,r_j-1),4)
        D4_eta(r_i,r_j,4) = -Dnonlinear_plus + Dnonlinear_minus
    end if
    !
    !
    call deallocate_vars_non_linear
    !
    !
end if
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
        !
        !
end subroutine compute_residue