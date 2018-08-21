! boundary_conditions_curv.f90
!
! implementation of boundary conditions
!
subroutine boundary_conditions_curv
    use vars
    implicit none
!
! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K
!
! inlet boundary
! 
    a_cr = (2.0d0*gama)*((gama-1.0d0)/(gama+1.0d0))*c_v*T_total
!
!
i = 1
do j = 2, jmax - 1
    T(i,j)   = T_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*(1.0d0+dtan(theta)**2)*((u(i,j)/a_cr)**2.0d0))
    p(i,j)   = p_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*(1.0d0+dtan(theta)**2)*((u(i,j)/a_cr)**2.0d0))**(gama/(gama-1.0d0))
    !
    ! para calcular a velocidade vamos usar as propriedades de estagnacao
    !
    Q_barra(i,j,1) = metric_jacobian(i,j)*(p(i,j)/(R*T(i,j)))
    a(i,j)   = sqrt(gama*p(i,j)*metric_jacobian(i,j)/Q_barra(i,j,1))
    u(i,j)   = sqrt( (2.0d0/(gama - 1.0d0))*( (T_total/T(i,j)) - 1.0d0 )*a(i,j)**2.0d0 )
    v(i,j)   = u(i,j)*dtan(theta)   
    Q_barra(i,j,2) = Q_barra(i,j,1)*u(i,j)
    Q_barra(i,j,3) = Q_barra(i,j,1)*v(i,j)    
    Q_barra(i,j,4) = p(i,j)/(gama-1.0d0) + (Q_barra(i,j,1)/(2.0d0*metric_jacobian(i,j)))*(u(i,j)**2 + v(i,j)**2)
    U_contravariant(i,j) = u(i,j)*ksi_x(i,j) + v(i,j)*ksi_y(i,j)
    V_contravariant(i,j) = u(i,j)*eta_x(i,j) + v(i,j)*eta_y(i,j)
end do
!
! outlet boundary
!
i = imax
do j = 2, jmax - 1 
    a(i,j) = sqrt(gama*p(i,j)*metric_jacobian(i,j)/Q_barra(i,j,1))
    q_vel(i,j) = sqrt(u(i,j)**2.0d0 + v(i,j)**2.0d0)
    if( (q_vel(i,j)/a(i,j)) < 1.0d0 ) then
        Q_barra(i,j,1) = Q_barra(i-1,j,1)
        Q_barra(i,j,2) = Q_barra(i-1,j,2)
        Q_barra(i,j,3) = Q_barra(i-1,j,3)
        p(i,j)   = p_total/3.0d0
        Q_barra(i,j,4) = p(i,j)/(gama-1.0d0) + (Q_barra(i,j,1)/(2.0d0*metric_jacobian(i,j)))*(u(i,j)**2 + v(i,j)**2)
    else
        Q_barra(i,j,1) = Q_barra(i-1,j,1)
        Q_barra(i,j,2) = Q_barra(i-1,j,2)
        Q_barra(i,j,3) = Q_barra(i-1,j,3)
        Q_barra(i,j,4) = Q_barra(i-1,j,4)
    end if
end do
!
! lower boundary 
!
j = 1
do i = 2, imax - 1
        U_contravariant(i,j) = U_contravariant(i,j+1)
        V_contravariant(i,j) = 0.0d0
        !
        ! extrapola os termos de metrica - duvida
        ! tem uma raiz dos termos de metrica - duvida
        ! 
        ! u e v sao determinados pela matriz jacobiana de transformacao
        ! so lebrar dos termos contrvariante das variaveis
        !
        u(i,j)               = x_ksi(i,j)*U_contravariant(i,j) 
        v(i,j)               = y_ksi(i,j)*U_contravariant(i,j) 
end do
!
!
!
j = 1
do i = 2, imax - 1
    Q_barra(i,j,1) = Q_barra(i,j+1,1)
    Q_barra(i,j,2) = Q_barra(i,j+1,2)
    Q_barra(i,j,3) = Q_barra(i,j+1,3)
    Q_barra(i,j,4) = p(i,j)/(gama-1.0d0) + (Q_barra(i,j,1)/2.0d0)*(u(i,j)**2 + v(i,j)**2) 
end do
!
! symmetry boundary
!
j = jmax
do i = 2, imax - 1
    Q_barra(i,j,1) = Q_barra(i,j-2,1)
    Q_barra(i,j,2) = Q_barra(i,j-2,2)
    Q_barra(i,j,3) =-Q_barra(i,j-2,3)
    Q_barra(i,j,4) = Q_barra(i,j-2,4)
end do
!
!
!
end subroutine boundary_conditions_curv
