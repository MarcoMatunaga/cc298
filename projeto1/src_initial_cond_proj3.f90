subroutine initial_cond_proj3
    use vars
    use vars_proj3
    implicit none
    integer(4)                      :: ind_i, ind_j
    real(8)                         :: u, v, p, rho

    p = 1.0d0/gama
    rho = 1.0d0
    u = mach_inlet
    v = 0.0d0

    do ind_i = 1, imax
        do ind_j = 1, jmax
            Q_barra(ind_i,ind_j,1) = metric_jacobian(ind_i,ind_j)*rho
            Q_barra(ind_i,ind_j,2) = metric_jacobian(ind_i,ind_j)*rho*u
            Q_barra(ind_i,ind_j,3) = metric_jacobian(ind_i,ind_j)*rho*v
            Q_barra(ind_i,ind_j,4) = metric_jacobian(ind_i,ind_j)*(p/(gama-1.0d0) + 0.50d0*rho*(u**2.0d0 + v**2.0d0))
        end do
    end do   

end subroutine initial_cond_proj3