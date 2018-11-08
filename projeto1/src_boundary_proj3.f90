subroutine boundary_proj3
    use vars
    use vars_proj3
    implicit none
    real(8)                         :: a_p3
    integer(4)                      :: bc_i, bc_j
    ! inlet boundary 
    bc_i = 1

    do bc_j = 1, jmax

        u_p3 = mach_inlet
        v_p3 = 0.0d0
        p_p3 = 1.0d0/gama
        rho  = 1.0d0

        Q_barra(bc_i,bc_j,1) = metric_jacobian(bc_i,bc_j)*rho
        Q_barra(bc_i,bc_j,2) = metric_jacobian(bc_i,bc_j)*rho*u_p3
        Q_barra(bc_i,bc_j,3) = metric_jacobian(bc_i,bc_j)*rho*v_p3
        Q_barra(bc_i,bc_j,4) = metric_jacobian(bc_i,bc_j)*(p_p3/(gama-1.0d0) &
                               + 0.5d0*rho*(u_p3**2.0d0 + v_p3**2.0d0))

    end do

    ! wall boundary
    bc_j = 1

        do bc_i = 2, imax - 1

        v_p3 = 0.0d0
        u_p3 = Q_barra(bc_i,bc_j+1,2)/Q_barra(bc_i,bc_j+1,1)
        rho  = Q_barra(bc_i,bc_j+1,1)/metric_jacobian(bc_i,bc_j)
        p_p3 = (gama-1.0d0) * (Q_barra(bc_i,bc_j+1,4)/metric_jacobian(bc_i,bc_j) & 
               - 0.5d0*( (Q_barra(bc_i,bc_j+1,2)/metric_jacobian(bc_i,bc_j))**2.0d0 &
               + (Q_barra(bc_i,bc_j+1,3)/metric_jacobian(bc_i,bc_j))**2.0d0)/rho)

        Q_barra(bc_i,bc_j,1) = metric_jacobian(bc_i,bc_j)*rho
        Q_barra(bc_i,bc_j,2) = metric_jacobian(bc_i,bc_j)*rho*u_p3
        Q_barra(bc_i,bc_j,3) = metric_jacobian(bc_i,bc_j)*rho*v_p3
        Q_barra(bc_i,bc_j,4) = metric_jacobian(bc_i,bc_j)*(p_p3/(gama-1.0d0) & 
                               + 0.5d0*rho*(u_p3**2.0d0 + v_p3**2.0d0))

    end do

    ! post shock conditions - upper boundary
        rho  =  (gama + 1.0d0)*mach_inlet_n**2.0d0/(2.0d0 + (gama - 1.0d0)*mach_inlet_n**2.0d0)
        p_p3 =  (1.0d0 + (2.0d0*gama/( gama + 1.0d0 ))*(mach_inlet_n**2.0d0 - 1.0d0) ) / gama
        a_p3 =  (gama*p_p3/rho)**0.5d0
        u_p3 =  mach_2 * a_p3 * cos(turn_angle)
        v_p3 = -mach_2 * a_p3 * sin(turn_angle)
    
    bc_j = jmax
    do bc_i = 2, imax - 1
        
        Q_barra(bc_i,bc_j,1) = metric_jacobian(bc_i,bc_j)*rho
        Q_barra(bc_i,bc_j,2) = metric_jacobian(bc_i,bc_j)*rho*u_p3
        Q_barra(bc_i,bc_j,3) = metric_jacobian(bc_i,bc_j)*rho*v_p3
        Q_barra(bc_i,bc_j,4) = metric_jacobian(bc_i,bc_j)*(p_p3/(gama-1.0d0) &
                               + 0.5d0*rho*(u_p3**2.0d0 + v_p3**2.0d0))

    end do

    ! outlet boundary
    bc_i = imax
    do bc_j = 1, jmax

        Q_barra(bc_i,bc_j,1) = Q_barra(bc_i-1,bc_j,1)
        Q_barra(bc_i,bc_j,2) = Q_barra(bc_i-1,bc_j,2)
        Q_barra(bc_i,bc_j,3) = Q_barra(bc_i-1,bc_j,3)
        Q_barra(bc_i,bc_j,4) = Q_barra(bc_i-1,bc_j,4)

    end do

end subroutine boundary_proj3