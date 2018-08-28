subroutine fluxes_curvilinear
    use vars
    implicit none
!
!  fluxes in curvilinear coordinates - ksi direction
!
do j = 1, jmax 
    do i = 1, imax 
        T(i,j)= T_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*(1.0d0+dtan(theta)**2)*((u(i,j)/a_cr)**2.0d0))
        p(i,j)= p_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*(1.0d0+dtan(theta)**2)*((u(i,j)/a_cr)**2.0d0))**(gama/(gama-1.0d0))
        Q_barra(i,j,1) = metric_jacobian(i,j)*(p(i,j)/(R*T(i,j)))
        u(i,j)= Q_barra(i,j,2)/Q_barra(i,j,1)   
        v(i,j)= Q_barra(i,j,3)/Q_barra(i,j,1) 
        Q_barra(i,j,4) = p(i,j)/(gama-1.0d0) + (Q_barra(i,j,1)/(2.0d0*metric_jacobian(i,j)))*(u(i,j)**2.0d0 + v(i,j)**2.0d0)
        U_contravariant(i,j) = u(i,j)*ksi_x(i,j) + v(i,j)*ksi_y(i,j)
        V_contravariant(i,j) = u(i,j)*eta_x(i,j) + v(i,j)*eta_y(i,j)
    end do
end do
!
!
j = jmax - 1
do i = i, imax
        V_contravariant(i,j) = 0.0d0
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
