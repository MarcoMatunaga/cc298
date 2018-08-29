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
    !Q_barra(i,j,4) = Q_barra(i,j,1)*( c_v*T(i,j) + 0.5d0*(u(i,j)**2.0d0 + v(i,j)**2.0d0) )
end do
!
! outlet boundary
!
i = imax
do j = 2, jmax - 1
    u(i,j)     = Q_barra(i-1,j,2)/Q_barra(i-1,j,1)
    v(i,j)     = Q_barra(i-1,j,3)/Q_barra(i-1,j,1)
    a(i,j)     = sqrt(gama*p(i-1,j)*metric_jacobian(i-1,j)/Q_barra(i-1,j,1))
    q_vel(i,j) = sqrt(u(i,j)**2.0d0 + v(i,j)**2.0d0)
    if( (q_vel(i,j)/a(i,j)) < 1.0d0 ) then
        write(*,*) "sub" , q_vel(i,j)
        Q_barra(i,j,1) = metric_jacobian(i,j)*(Q_barra(i-1,j,1)/metric_jacobian(i-1,j))
        Q_barra(i,j,2) = metric_jacobian(i,j)*(Q_barra(i-1,j,2)/metric_jacobian(i-1,j))
        Q_barra(i,j,3) = metric_jacobian(i,j)*(Q_barra(i-1,j,3)/metric_jacobian(i-1,j))
        p(i,j)         = p_total/3.0d0
        Q_barra(i,j,4) = metric_jacobian(i,j)*(p(i,j)/(gama-1.0d0) &
                        + (Q_barra(i,j,1)/(2.0d0*metric_jacobian(i,j)))*(u(i,j)**2.0d0 + v(i,j)**2.0d0))
    else
        write(*,*) "super", q_vel(i,j)
        Q_barra(i,j,1) = metric_jacobian(i,j)*(Q_barra(i-1,j,1)/metric_jacobian(i-1,j))
        Q_barra(i,j,2) = metric_jacobian(i,j)*(Q_barra(i-1,j,2)/metric_jacobian(i-1,j))
        Q_barra(i,j,3) = metric_jacobian(i,j)*(Q_barra(i-1,j,3)/metric_jacobian(i-1,j))
        Q_barra(i,j,4) = metric_jacobian(i,j)*(Q_barra(i-1,j,4)/metric_jacobian(i-1,j))
    end if
end do
!
! lower boundary (wall boundary)
!
j = 1
do i = 1, imax
    U_contravariant(i,j) = U_contravariant(i,j+1)
    V_contravariant(i,j) = 0.0d0
        !
        ! extrapola os termos de metrica - duvida
        ! tem uma raiz dos termos de metrica - duvida
        ! 
        ! u e v sao determinados pela matriz jacobiana de transformacao
        ! so lebrar dos termos contrvariante das variaveis
        !
    u(i,j)         = x_ksi(i,j)*U_contravariant(i,j) 
    v(i,j)         = y_ksi(i,j)*U_contravariant(i,j) 
    Q_barra(i,j,1) = metric_jacobian(i,j)*(Q_barra(i,j+1,1)/metric_jacobian(i,j+1))
    Q_barra(i,j,2) = metric_jacobian(i,j)*(Q_barra(i,j+1,2)/metric_jacobian(i,j+1))
    Q_barra(i,j,3) = metric_jacobian(i,j)*(Q_barra(i,j+1,3)/metric_jacobian(i,j+1))
    Q_barra(i,j,4) = metric_jacobian(i,j)*(Q_barra(i,j+1,4)/metric_jacobian(i,j+1))
end do
!
! symmetry boundary
!
j = jmax
do i = 1, imax 
    Q_barra(i,j,1) =  metric_jacobian(i,j)*(Q_barra(i,j-2,1)/metric_jacobian(i,j-2))
    Q_barra(i,j,2) =  metric_jacobian(i,j)*(Q_barra(i,j-2,2)/metric_jacobian(i,j-2))
    Q_barra(i,j,3) = -metric_jacobian(i,j)*(Q_barra(i,j-2,3)/metric_jacobian(i,j-2))
    Q_barra(i,j,4) =  metric_jacobian(i,j)*(Q_barra(i,j-2,4)/metric_jacobian(i,j-2))
end do
!
!
!
! do j = 1, jmax
!     do i = 1, imax
!         if (p(i,j) < 0.00003d0) then 
!                 print *, i,j,iter,p(i,j)
!         end if
!     end do     
! end do
!
!
!
end subroutine boundary_conditions_curv
