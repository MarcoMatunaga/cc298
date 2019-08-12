subroutine artificial_dissipation(i_diss,j_diss) 
    use vars
    use functions
    implicit none
integer(4)                      :: i_diss, j_diss
!
!
if (which_diss == 1) then 
    !
    !
    D4_ksi(i_diss,j_diss,1) = dis_ksi4(i_diss,j_diss,1,eps_dis_e)
    D4_ksi(i_diss,j_diss,2) = dis_ksi4(i_diss,j_diss,2,eps_dis_e)
    D4_ksi(i_diss,j_diss,3) = dis_ksi4(i_diss,j_diss,3,eps_dis_e)
    D4_ksi(i_diss,j_diss,4) = dis_ksi4(i_diss,j_diss,4,eps_dis_e)
    !
    !
    D4_eta(i_diss,j_diss,1) = dis_eta4(i_diss,j_diss,1,eps_dis_e)
    D4_eta(i_diss,j_diss,2) = dis_eta4(i_diss,j_diss,2,eps_dis_e)
    D4_eta(i_diss,j_diss,3) = dis_eta4(i_diss,j_diss,3,eps_dis_e)
    D4_eta(i_diss,j_diss,4) = dis_eta4(i_diss,j_diss,4,eps_dis_e)
    !
    !
else if (which_diss == 2) then 
    !
    !
    D4_ksi(i_diss,j_diss,1) = dis_ksi2(i_diss,j_diss,1,eps_dis_e)
    D4_ksi(i_diss,j_diss,2) = dis_ksi2(i_diss,j_diss,2,eps_dis_e)
    D4_ksi(i_diss,j_diss,3) = dis_ksi2(i_diss,j_diss,3,eps_dis_e)
    D4_ksi(i_diss,j_diss,4) = dis_ksi2(i_diss,j_diss,4,eps_dis_e)
    !
    !
    D4_eta(i_diss,j_diss,1) = dis_eta2(i_diss,j_diss,1,eps_dis_e)
    D4_eta(i_diss,j_diss,2) = dis_eta2(i_diss,j_diss,2,eps_dis_e)
    D4_eta(i_diss,j_diss,3) = dis_eta2(i_diss,j_diss,3,eps_dis_e)
    D4_eta(i_diss,j_diss,4) = dis_eta2(i_diss,j_diss,4,eps_dis_e)
    !
    !
else if (which_diss == 3) then 
    !
    !
    call allocate_vars_non_linear
    call non_linear_dissipation
    !
    !
    if (i_diss == 2 .or. i_diss == imax - 1) then
            D4_ksi(i_diss,j_diss,1) = dis_ksi2(i_diss,j_diss,1,eps_dis_e)
            D4_ksi(i_diss,j_diss,2) = dis_ksi2(i_diss,j_diss,2,eps_dis_e)
            D4_ksi(i_diss,j_diss,3) = dis_ksi2(i_diss,j_diss,3,eps_dis_e)
            D4_ksi(i_diss,j_diss,4) = dis_ksi2(i_diss,j_diss,4,eps_dis_e)
    else if (j_diss == 2 .or. j_diss == jmax - 1) then
            D4_eta(i_diss,j_diss,1) = dis_eta2(i_diss,j_diss,1,eps_dis_e)
            D4_eta(i_diss,j_diss,2) = dis_eta2(i_diss,j_diss,2,eps_dis_e)
            D4_eta(i_diss,j_diss,3) = dis_eta2(i_diss,j_diss,3,eps_dis_e)
            D4_eta(i_diss,j_diss,4) = dis_eta2(i_diss,j_diss,4,eps_dis_e)
    else 
        Dnonlinear_plus   = non_lin_dis_ksi(i_diss,j_diss,eps2_ksi(i_diss,j_diss),eps4_ksi(i_diss,j_diss),1)
        Dnonlinear_minus  = non_lin_dis_ksi(i_diss-1,j_diss,eps2_ksi(i_diss-1,j_diss),eps4_ksi(i_diss-1,j_diss),1)
        D4_ksi(i_diss,j_diss,1) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_ksi(i_diss,j_diss,eps2_ksi(i_diss,j_diss),eps4_ksi(i_diss,j_diss),2)
        Dnonlinear_minus  = non_lin_dis_ksi(i_diss-1,j_diss,eps2_ksi(i_diss-1,j_diss),eps4_ksi(i_diss-1,j_diss),2)
        D4_ksi(i_diss,j_diss,2) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_ksi(i_diss,j_diss,eps2_ksi(i_diss,j_diss),eps4_ksi(i_diss,j_diss),3)
        Dnonlinear_minus  = non_lin_dis_ksi(i_diss-1,j_diss,eps2_ksi(i_diss-1,j_diss),eps4_ksi(i_diss-1,j_diss),3)
        D4_ksi(i_diss,j_diss,3) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_ksi(i_diss,j_diss,eps2_ksi(i_diss,j_diss),eps4_ksi(i_diss,j_diss),4)
        Dnonlinear_minus  = non_lin_dis_ksi(i_diss-1,j_diss,eps2_ksi(i_diss-1,j_diss),eps4_ksi(i_diss-1,j_diss),4)
        D4_ksi(i_diss,j_diss,4) = -Dnonlinear_plus + Dnonlinear_minus
        !
        !
        Dnonlinear_plus   = non_lin_dis_eta(i_diss,j_diss,eps2_eta(i_diss,j_diss),eps4_eta(i_diss,j_diss),1)
        Dnonlinear_minus  = non_lin_dis_eta(i_diss,j_diss-1,eps2_eta(i_diss,j_diss-1),eps4_eta(i_diss,j_diss-1),1)
        D4_eta(i_diss,j_diss,1) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_eta(i_diss,j_diss,eps2_eta(i_diss,j_diss),eps4_eta(i_diss,j_diss),2)
        Dnonlinear_minus  = non_lin_dis_eta(i_diss,j_diss-1,eps2_eta(i_diss,j_diss-1),eps4_eta(i_diss,j_diss-1),2)
        D4_eta(i_diss,j_diss,2) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_eta(i_diss,j_diss,eps2_eta(i_diss,j_diss),eps4_eta(i_diss,j_diss),3)
        Dnonlinear_minus  = non_lin_dis_eta(i_diss,j_diss-1,eps2_eta(i_diss,j_diss-1),eps4_eta(i_diss,j_diss-1),3)
        D4_eta(i_diss,j_diss,3) = -Dnonlinear_plus + Dnonlinear_minus
        !
        Dnonlinear_plus   = non_lin_dis_eta(i_diss,j_diss,eps2_eta(i_diss,j_diss),eps4_eta(i_diss,j_diss),4)
        Dnonlinear_minus  = non_lin_dis_eta(i_diss,j_diss-1,eps2_eta(i_diss,j_diss-1),eps4_eta(i_diss,j_diss-1),4)
        D4_eta(i_diss,j_diss,4) = -Dnonlinear_plus + Dnonlinear_minus
    end if
    !
    !
    call deallocate_vars_non_linear
    !
    !
end if
!
!
end subroutine artificial_dissipation