
! implementation of boundary conditions

subroutine boundary_conditions_curv
    use vars
    implicit none

! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K

! inlet boundary

i = 1
do j = 1, jmax 
    u(i,j)  = Q_barra(i+1,j,2)/Q_barra(i+1,j,1)
    v(i,j)  = u(i,j)*dtan(theta)  
    T(i,j)  = T_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*((u(i,j)/a_cr)**2.0d0))
    p(i,j)  = p_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*((u(i,j)/a_cr)**2.0d0))**(gama/(gama-1.0d0))
    Q_barra(i,j,1) = metric_jacobian(i,j)*(p(i,j)/(R*T(i,j)))
    Q_barra(i,j,2) = Q_barra(i,j,1)*u(i,j)
    Q_barra(i,j,3) = Q_barra(i,j,1)*v(i,j)    
    Q_barra(i,j,4) = metric_jacobian(i,j)*(p(i,j)/(gama-1.0d0) &
                     + (Q_barra(i,j,1)/(2.0d0*metric_jacobian(i,j)))*(u(i,j)**2.0d0 + v(i,j)**2.0d0))

    ! E_barra(i,j,1) = Q_barra(i,j,1)*U_contravariant(i,j) 
    ! E_barra(i,j,2) = Q_barra(i,j,2)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_x(i,j)
    ! E_barra(i,j,3) = Q_barra(i,j,3)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_y(i,j)
    ! E_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*U_contravariant(i,j)

    ! F_barra(i,j,1) = Q_barra(i,j,1)*V_contravariant(i,j) 
    ! F_barra(i,j,2) = Q_barra(i,j,2)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_x(i,j)
    ! F_barra(i,j,3) = Q_barra(i,j,3)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_y(i,j)
    ! F_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*V_contravariant(i,j)
end do

! lower boundary (wall boundary)

j = 1
do i = 2, imax
    U_contravariant(i,j) = U_contravariant(i,j+1)
    V_contravariant(i,j) = 0.0d0
        
        ! u e v sao determinados pela matriz jacobiana de transformacao
        ! so lebrar dos termos contrvariante das variaveis

    u(i,j)         = x_ksi(i,j)*U_contravariant(i,j) 
    v(i,j)         = y_ksi(i,j)*U_contravariant(i,j) 
    Q_barra(i,j,1) = metric_jacobian(i,j)*(Q_barra(i,j+1,1)/metric_jacobian(i,j+1))
    Q_barra(i,j,2) = metric_jacobian(i,j)*(Q_barra(i,j+1,2)/metric_jacobian(i,j+1))
    Q_barra(i,j,3) = metric_jacobian(i,j)*(Q_barra(i,j+1,3)/metric_jacobian(i,j+1))
    Q_barra(i,j,4) = metric_jacobian(i,j)*(Q_barra(i,j+1,4)/metric_jacobian(i,j+1))

    ! E_barra(i,j,1) = Q_barra(i,j,1)*U_contravariant(i,j) 
    ! E_barra(i,j,2) = Q_barra(i,j,2)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_x(i,j)
    ! E_barra(i,j,3) = Q_barra(i,j,3)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_y(i,j)
    ! E_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*U_contravariant(i,j)

    ! F_barra(i,j,1) = Q_barra(i,j,1)*V_contravariant(i,j) 
    ! F_barra(i,j,2) = Q_barra(i,j,2)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_x(i,j)
    ! F_barra(i,j,3) = Q_barra(i,j,3)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_y(i,j)
    ! F_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*V_contravariant(i,j)
end do

! symmetry boundary

j = jmax
do i = 1, imax
    Q_barra(i,j,1) =  metric_jacobian(i,j)*(Q_barra(i,j-2,1)/metric_jacobian(i,j-2))
    Q_barra(i,j,2) =  metric_jacobian(i,j)*(Q_barra(i,j-2,2)/metric_jacobian(i,j-2))
    Q_barra(i,j,3) = -metric_jacobian(i,j)*(Q_barra(i,j-2,3)/metric_jacobian(i,j-2))
    Q_barra(i,j,4) =  metric_jacobian(i,j)*(Q_barra(i,j-2,4)/metric_jacobian(i,j-2))

    ! E_barra(i,j,1) = Q_barra(i,j,1)*U_contravariant(i,j) 
    ! E_barra(i,j,2) = Q_barra(i,j,2)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_x(i,j)
    ! E_barra(i,j,3) = Q_barra(i,j,3)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_y(i,j)
    ! E_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*U_contravariant(i,j)

    ! F_barra(i,j,1) = Q_barra(i,j,1)*V_contravariant(i,j) 
    ! F_barra(i,j,2) = Q_barra(i,j,2)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_x(i,j)
    ! F_barra(i,j,3) = Q_barra(i,j,3)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_y(i,j)
    ! F_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*V_contravariant(i,j)
end do

! outlet boundary

i = imax
do j = 1, jmax
    u(i,j)     = Q_barra(i,j,2)/Q_barra(i,j,1)
    v(i,j)     = Q_barra(i,j,3)/Q_barra(i,j,1)
    a(i,j)     = sqrt(gama*p(i,j)*metric_jacobian(i,j)/Q_barra(i,j,1))
    q_vel(i,j) = sqrt(u(i,j)**2.0d0 + v(i,j)**2.0d0)
    if( (q_vel(i,j)/a(i,j)) < 1.0d0 ) then
        Q_barra(i,j,1) = metric_jacobian(i,j)*(Q_barra(i-1,j,1)/metric_jacobian(i-1,j))
        Q_barra(i,j,2) = metric_jacobian(i,j)*(Q_barra(i-1,j,2)/metric_jacobian(i-1,j))
        Q_barra(i,j,3) = metric_jacobian(i,j)*(Q_barra(i-1,j,3)/metric_jacobian(i-1,j))
        p(i,j)         = p_total/3.0d0
        Q_barra(i,j,4) = metric_jacobian(i,j)*(p(i,j)/(gama-1.0d0) &
                        + (Q_barra(i,j,1)/(2.0d0*metric_jacobian(i,j)))*(u(i,j)**2.0d0 + v(i,j)**2.0d0))

        ! E_barra(i,j,1) = Q_barra(i,j,1)*U_contravariant(i,j) 
        ! E_barra(i,j,2) = Q_barra(i,j,2)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_x(i,j)
        ! E_barra(i,j,3) = Q_barra(i,j,3)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_y(i,j)
        ! E_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*U_contravariant(i,j)

        ! F_barra(i,j,1) = Q_barra(i,j,1)*V_contravariant(i,j) 
        ! F_barra(i,j,2) = Q_barra(i,j,2)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_x(i,j)
        ! F_barra(i,j,3) = Q_barra(i,j,3)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_y(i,j)
        ! F_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*V_contravariant(i,j)
    else
        Q_barra(i,j,1) = metric_jacobian(i,j)*(Q_barra(i-1,j,1)/metric_jacobian(i-1,j))
        Q_barra(i,j,2) = metric_jacobian(i,j)*(Q_barra(i-1,j,2)/metric_jacobian(i-1,j))
        Q_barra(i,j,3) = metric_jacobian(i,j)*(Q_barra(i-1,j,3)/metric_jacobian(i-1,j))
        Q_barra(i,j,4) = metric_jacobian(i,j)*(Q_barra(i-1,j,4)/metric_jacobian(i-1,j))

        ! E_barra(i,j,1) = Q_barra(i,j,1)*U_contravariant(i,j) 
        ! E_barra(i,j,2) = Q_barra(i,j,2)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_x(i,j)
        ! E_barra(i,j,3) = Q_barra(i,j,3)*U_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*ksi_y(i,j)
        ! E_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*U_contravariant(i,j)

        ! F_barra(i,j,1) = Q_barra(i,j,1)*V_contravariant(i,j) 
        ! F_barra(i,j,2) = Q_barra(i,j,2)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_x(i,j)
        ! F_barra(i,j,3) = Q_barra(i,j,3)*V_contravariant(i,j) + metric_jacobian(i,j)*p(i,j)*eta_y(i,j)
        ! F_barra(i,j,4) = metric_jacobian(i,j)*( (Q_barra(i,j,4)/metric_jacobian(i,j)) + p(i,j) )*V_contravariant(i,j)
    end if
end do

! do j = 1, jmax
!     do i = 1, imax
!         if (p(i,j) < 0.00003d0) then 
!                 print *, i,j,iter,p(i,j)
!         end if
!     end do     
! end do

end subroutine boundary_conditions_curv
