subroutine non_linear_dissipation
    use vars
    implicit none
    !
    integer(4)                        :: i_nld, j_nld
    real(8)                           :: max_nu
    !
    !*****
    ! we can create teo new variables one
    ! each one for the jacobian matrix spectral radius 
    !*****
    !
    !******
    ! we can calculate for specific points we do not 
    ! use all points of the arrays below
    !******
    !
    ! the speed of sound came from another subroutine
    !
    !
    ! i_nld = 1
    ! do j_nld = 1, jmax
    !     nu_dis_ksi(i_nld,j_nld) = abs(p(i_nld,j_nld) - 2.0d0*p(i_nld,j_nld) + p(i_nld,j_nld))/(p(i_nld,j_nld) &
    !                               - 2.0d0*p(i_nld,j_nld) + p(i_nld,j_nld))
    ! end do 
    ! !
    ! !
    ! i_nld = imax
    ! do j_nld = 1, jmax
    !     nu_dis_ksi(i_nld,j_nld) =
    ! end do
    ! !
    ! !
    ! j = 1
    ! do i_nld = 1, imax 
    !     nu_dis_eta(i_nld,j_nld) =
    ! end do
    ! !
    ! !
    ! j = jmax
    ! do i_nld = 1, imax
    !     nu_dis_eta(i_nld,j_nld) =
    ! end do
    nu_dis_eta = 0.0d0
    nu_dis_ksi = 0.0d0
    !
    !
    do j_nld = 2, jmax - 1
        do i_nld = 2, imax - 1
            nu_dis_ksi(i_nld,j_nld) = abs(p(i_nld+1,j_nld) - 2.0d0*p(i_nld,j_nld) + p(i_nld-1,j_nld))/(p(i_nld+1,j_nld) &
                                      + 2.0d0*p(i_nld,j_nld) + p(i_nld-1,j_nld))     
            !
            !
            nu_dis_eta(i_nld,j_nld) = abs(p(i_nld,j_nld+1) - 2.0d0*p(i_nld,j_nld) + p(i_nld,j_nld-1))/(p(i_nld,j_nld+1) &
                                      + 2.0d0*p(i_nld,j_nld) + p(i_nld,j_nld-1))       
        end do
    end do
    !
    !
    do j_nld = 2, jmax - 1
        do i_nld = 2, imax - 1
            !
            !
            max_nu = max(nu_dis_ksi(i_nld+1,j_nld),nu_dis_ksi(i_nld,j_nld),nu_dis_ksi(i_nld-1,j_nld))
            eps2_ksi(i_nld,j_nld)   = k2*delta_t(i_nld,j_nld)*max_nu
            eps4_ksi(i_nld,j_nld)   = max(0.0d0,(k4*delta_t(i_nld,j_nld)-eps2_ksi(i_nld,j_nld)))
            !
            !
            max_nu = max(nu_dis_eta(i_nld,j_nld+1),nu_dis_eta(i_nld,j_nld),nu_dis_eta(i_nld,j_nld-1))
            eps2_eta(i_nld,j_nld)   = k2*delta_t(i_nld,j_nld)*max_nu
            eps4_eta(i_nld,j_nld)   = max(0.0d0,(k4*delta_t(i_nld,j_nld)-eps2_eta(i_nld,j_nld)))
            !
            !
        end do 
    end do
    !
    !
    do j_nld = 1, jmax
        do i_nld = 1, imax
                sigma_ksi(i_nld,j_nld)      = abs(U_contravariant(i_nld,j_nld)) & 
                                              + a(i_nld,j_nld)*sqrt(ksi_x(i_nld,j_nld)**2.0d0 + ksi_y(i_nld,j_nld)**2.0d0 )
                sigma_eta(i_nld,j_nld)      = abs(V_contravariant(i_nld,j_nld)) &
                                              + a(i_nld,j_nld)*sqrt(eta_x(i_nld,j_nld)**2.0d0 + eta_y(i_nld,j_nld)**2.0d0 )
        end do 
    end do
    !
    !******
    ! there are problems
    !       which direction use
    !       maximum or minimum value for dissipation
    !******
    !
    eps_dis_i  = dis_factor*(maxval(eps2_ksi)+maxvaL(eps4_ksi)) 
    ! write(*,*) eps_dis_i
    !
    !
end subroutine non_linear_dissipation