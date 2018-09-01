subroutine fluxes_curvilinear
    use vars
    implicit none
!
!
!
E_barra = 0.0d0
F_barra = 0.0d0
!
!  fluxes in curvilinear coordinates - ksi direction
!
do j = 1, jmax 
    do i = 1, imax
        u(i,j) = Q_barra(i,j,2)/Q_barra(i,j,1)   
        v(i,j) = Q_barra(i,j,3)/Q_barra(i,j,1)
        p(i,j) = (gama - 1.0d0)*( Q_barra(i,j,4)/metric_jacobian(i,j) &
                 - 0.50d0*( (Q_barra(i,j,1)/metric_jacobian(i,j))*(u(i,j)**2.0d0+v(i,j)**2.0d0) ) )
        U_contravariant(i,j) = u(i,j)*ksi_x(i,j) + v(i,j)*ksi_y(i,j)
        V_contravariant(i,j) = u(i,j)*eta_x(i,j) + v(i,j)*eta_y(i,j)
    if (j == 1) V_contravariant(i,j) = 0.0d0
    end do
end do
!
!
!
do j = 1, jmax
    do i = 1, imax
        E_barra(i,j,1) = Q_barra(i,j,1)*U_contravariant(i,j) 
        E_barra(i,j,2) = Q_barra(i,j,2)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_x(i,j)
        E_barra(i,j,3) = Q_barra(i,j,3)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_y(i,j)
        E_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*U_contravariant(i,j)
    end do 
end do
!
!  fluxes in curvilinear coordinates - eta direction
!
do j = 1, jmax
    do i = 1, imax
        F_barra(i,j,1) = Q_barra(i,j,1)*V_contravariant(i,j) 
        F_barra(i,j,2) = Q_barra(i,j,2)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_x(i,j)
        F_barra(i,j,3) = Q_barra(i,j,3)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_y(i,j)
        F_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*V_contravariant(i,j)
    end do         
end do 
!
!
! 
end subroutine fluxes_curvilinear
