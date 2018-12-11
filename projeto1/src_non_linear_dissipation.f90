subroutine non_linear_dissipation
    use vars
    implicit none
    !
    integer(4)                          :: i_nld, j_nld
    real(8)                             :: max_nu, a, u, v, rho
    real(8),dimension(:,:),allocatable  :: p_v
    
    allocate(p_v(imax,jmax))     
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
    !     nu_dis_ksi(i_nld,j_nld) = abs(p_v(i_nld,j_nld) - 2.0d0*p_v(i_nld,j_nld) + p_v(i_nld,j_nld))/(p_v(i_nld,j_nld) &
    !                               - 2.0d0*p_v(i_nld,j_nld) + p_v(i_nld,j_nld))
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

    do j_nld = 1, jmax 
        do i_nld = 1, imax 
            u   = Q_barra(i_nld,j_nld,2)/Q_barra(i_nld,j_nld,1)
            v   = Q_barra(i_nld,j_nld,3)/Q_barra(i_nld,j_nld,1)
            rho = Q_barra(i_nld,j_nld,1)/metric_jacobian(i_nld,j_nld)
            p_v(i_nld,j_nld)   = (gama-1.0d0)*(Q_barra(i_nld,j_nld,4)/metric_jacobian(i_nld,j_nld) &
                                 - 0.50d0*rho*(u**2.0d0+v**2.0d0))
        end do 
    end do

    do j_nld = 2, jmax - 1
        do i_nld = 2, imax - 1
            nu_dis_ksi(i_nld,j_nld) = abs(p_v(i_nld+1,j_nld) - 2.0d0*p_v(i_nld,j_nld) + p_v(i_nld-1,j_nld))/(p_v(i_nld+1,j_nld) &
                                      + 2.0d0*p_v(i_nld,j_nld) + p_v(i_nld-1,j_nld))     

            nu_dis_eta(i_nld,j_nld) = abs(p_v(i_nld,j_nld+1) - 2.0d0*p_v(i_nld,j_nld) + p_v(i_nld,j_nld-1))/(p_v(i_nld,j_nld+1) &
                                      + 2.0d0*p_v(i_nld,j_nld) + p_v(i_nld,j_nld-1))       
        end do
    end do

    do j_nld = 2, jmax - 1
        do i_nld = 2, imax - 1

            max_nu = max(nu_dis_ksi(i_nld+1,j_nld),nu_dis_ksi(i_nld,j_nld),nu_dis_ksi(i_nld-1,j_nld))
            eps2_ksi(i_nld,j_nld)   = k2*delta_t(i_nld,j_nld)*max_nu
            eps4_ksi(i_nld,j_nld)   = max(0.0d0,(k4*delta_t(i_nld,j_nld)-eps2_ksi(i_nld,j_nld)))

            max_nu = max(nu_dis_eta(i_nld,j_nld+1),nu_dis_eta(i_nld,j_nld),nu_dis_eta(i_nld,j_nld-1))
            eps2_eta(i_nld,j_nld)   = k2*delta_t(i_nld,j_nld)*max_nu
            eps4_eta(i_nld,j_nld)   = max(0.0d0,(k4*delta_t(i_nld,j_nld)-eps2_eta(i_nld,j_nld)))

        end do 
    end do

    do j_nld = 1, jmax
        do i_nld = 1, imax
                a  = sqrt(gama*p_v(i_nld,j_nld)/(Q_barra(i_nld,j_nld,1)/metric_jacobian(i_nld,j_nld)))
                sigma_ksi(i_nld,j_nld)      = abs(U_contravariant(i_nld,j_nld)) & 
                                              + a*sqrt(ksi_x(i_nld,j_nld)**2.0d0 + ksi_y(i_nld,j_nld)**2.0d0 )
                sigma_eta(i_nld,j_nld)      = abs(V_contravariant(i_nld,j_nld)) &
                                              + a*sqrt(eta_x(i_nld,j_nld)**2.0d0 + eta_y(i_nld,j_nld)**2.0d0 )
        end do 
    end do
    !
    !******
    ! there are problems
    !       which direction use
    !       maximum or minimum value for dissipation
    !******
    !
    eps_dis_i = dis_factor*(min(minval(eps2_ksi),minval(eps2_eta))+min(minval(eps4_eta),minvaL(eps4_ksi))) 

    deallocate(p_v)
end subroutine non_linear_dissipation