subroutine initial_conditions_curv
    use vars
    implicit none
!
! inlet boundary
! 
! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K
!
T = T_total
p = p_total
u = 0.0d0
v = 0.0d0
do j = 1, jmax
    do i = 1, imax
            Q_barra(i,j,1) = 0.0d0
            Q_barra(i,j,2) = 0.0d0
            Q_barra(i,j,3) = 0.0d0
            Q_barra(i,j,4) = 0.0d0
    end do        
end do
!
! Euler vectors in cartesian coordinates
!
allocate(Q(imax,jmax,dim))
!
!
!
do j = 1, jmax
    do i = 1, imax
        Q(i,j,1) = p(i,j)/(R*T(i,j))    
        Q(i,j,2) = Q(i,j,1)*u(i,j)
        Q(i,j,3) = Q(i,j,1)*v(i,j)
        Q(i,j,4) = p(i,j)/(gama-1.0d0) + (Q(i,j,1)/2.0d0)*(u(i,j)**2 + v(i,j)**2)
    end do
end do
!
! 
!
do j = 1, jmax
    do i = 1, imax
        Q_barra(i,j,1) = metric_jacobian(i,j)*Q(i,j,1)
        Q_barra(i,j,2) = metric_jacobian(i,j)*Q(i,j,2)
        Q_barra(i,j,3) = metric_jacobian(i,j)*Q(i,j,3)
        Q_barra(i,j,4) = metric_jacobian(i,j)*Q(i,j,4)
    end do
end do  
!
! change the outlet values
!
i = imax
do j = 1, jmax
    p(i,j)         = p_total/3.0d0
    Q(i,j,1)       = p(i,j)/(R*T(i,j))
    Q(i,j,4)       = p(i,j)/(gama-1.0d0) + (Q(i,j,1)/2.0d0)*(u(i,j)**2 + v(i,j)**2)
    Q_barra(i,j,1) = metric_jacobian(i,j)*Q(i,j,1)
    Q_barra(i,j,4) = metric_jacobian(i,j)*Q(i,j,4)
end do
!
!  fluxes in curvilinear coordinates - ksi direction
!
do j = 1, jmax
    do i = 1, imax
        U_contravariant(i,j) = u(i,j)*ksi_x(i,j) + v(i,j)*ksi_y(i,j)
        E_barra(i,j,1) = Q_barra(i,j,1)
        E_barra(i,j,2) = Q_barra(i,j,2)*U_contravariant(i,j) + p(i,j)*ksi_x(i,j)
        E_barra(i,j,3) = Q_barra(i,j,3)*U_contravariant(i,j) + p(i,j)*ksi_y(i,j)
        E_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*( U_contravariant(i,j))
    end do 
end do
!
!  fluxes in curvilinear coordinates - eta direction
!
do j = 1, jmax
    do i = 1, imax
        V_contravariant(i,j) = u(i,j)*eta_x(i,j) + v(i,j)*eta_y(i,j)
        F_barra(i,j,1) = Q_barra(i,j,1)
        F_barra(i,j,2) = Q_barra(i,j,2)*V_contravariant(i,j) + p(i,j)*eta_x(i,j)
        F_barra(i,j,3) = Q_barra(i,j,3)*V_contravariant(i,j) + p(i,j)*eta_y(i,j)
        F_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*( V_contravariant(i,j))
    end do         
end do 
!
!
! 
deallocate(Q)
!
!
end subroutine initial_conditions_curv