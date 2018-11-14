subroutine fluxes_curvilinear
    use vars
    use fluxes_pos_neg
    implicit none
    real(8)                         :: u, v, p
!********
!modify: calculate only in the inner points 
! e pass the boundary points to the boundary subroutine
!********
!********
!change the way you calculate
!********
E_barra = 0.0d0
F_barra = 0.0d0

do j = 1, jmax
    do i = 1, imax
               
        u = Q_barra(i,j,2)/Q_barra(i,j,1)   
        v = Q_barra(i,j,3)/Q_barra(i,j,1)
        p = (gama - 1.0d0)*( Q_barra(i,j,4)/metric_jacobian(i,j) &
            - 0.50d0*( (Q_barra(i,j,1)/metric_jacobian(i,j))*(u**2.0d0+v**2.0d0) ))

        if ( j /= 1 .and. j/= jmax .and. i /= 1 .and. i/= imax ) then
            U_contravariant(i,j) = u*ksi_x(i,j) + v*ksi_y(i,j)
            V_contravariant(i,j) = u*eta_x(i,j) + v*eta_y(i,j)
        endif

        !  fluxes in curvilinear coordinates - ksi direction
        E_barra(i,j,1) = Q_barra(i,j,1)*U_contravariant(i,j) 
        E_barra(i,j,2) = Q_barra(i,j,2)*U_contravariant(i,j) + metric_jacobian(i,j)*p*ksi_x(i,j)
        E_barra(i,j,3) = Q_barra(i,j,3)*U_contravariant(i,j) + metric_jacobian(i,j)*p*ksi_y(i,j)
        E_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p )*U_contravariant(i,j)
        
        !  fluxes in curvilinear coordinates - eta direction
        F_barra(i,j,1) = Q_barra(i,j,1)*V_contravariant(i,j) 
        F_barra(i,j,2) = Q_barra(i,j,2)*V_contravariant(i,j) + metric_jacobian(i,j)*p*eta_x(i,j)
        F_barra(i,j,3) = Q_barra(i,j,3)*V_contravariant(i,j) + metric_jacobian(i,j)*p*eta_y(i,j)
        F_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p )*V_contravariant(i,j)

    end do 
end do

end subroutine fluxes_curvilinear