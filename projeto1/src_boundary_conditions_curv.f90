
! implementation of boundary conditions

subroutine boundary_conditions_curv
    use vars
    implicit none
        real(8)                         :: u, v, p, rho, T, a, q_vel

! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K

! inlet boundary

i = 1
do j = 2, jmax - 1
    u  = Q_barra(i+1,j,2)/Q_barra(i+1,j,1)
    v  = u*dtan(theta)  
    T  = T_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*((u/a_cr)**2.0d0))
    p  = p_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*((u/a_cr)**2.0d0))**(gama/(gama-1.0d0))
    Q_barra(i,j,1) = metric_jacobian(i,j)*(p/(R*T))
    Q_barra(i,j,2) = Q_barra(i,j,1)*u
    Q_barra(i,j,3) = Q_barra(i,j,1)*v    
    Q_barra(i,j,4) = metric_jacobian(i,j)*(p/(gama-1.0d0) &
                     + (Q_barra(i,j,1)/(2.0d0*metric_jacobian(i,j)))*(u**2.0d0 + v**2.0d0))
end do

! lower boundary (wall boundary)

j = 1
do i = 1, imax
    U_contravariant(i,j) = U_contravariant(i,j+1)
    V_contravariant(i,j) = 0.0d0
        
        ! u e v sao determinados pela matriz jacobiana de transformacao
        ! so lebrar dos termos contrvariante das variaveis
    Q_barra(i,j,1) = metric_jacobian(i,j)*(Q_barra(i,j+1,1)/metric_jacobian(i,j+1))
    Q_barra(i,j,2) = metric_jacobian(i,j)*(Q_barra(i,j+1,2)/metric_jacobian(i,j+1))
    Q_barra(i,j,3) = metric_jacobian(i,j)*(Q_barra(i,j+1,3)/metric_jacobian(i,j+1))
    Q_barra(i,j,4) = metric_jacobian(i,j)*(Q_barra(i,j+1,4)/metric_jacobian(i,j+1))
end do

! outlet boundary

i = imax
do j = 2, jmax - 1
    rho = Q_barra(i,j,1)/metric_jacobian(i,j)
    u   = Q_barra(i,j,2)/Q_barra(i,j,1)
    v   = Q_barra(i,j,3)/Q_barra(i,j,1)
    p = (gama-1.0d0) * (Q_barra(i,j,4)/metric_jacobian(i,j) & 
               - 0.5d0*rho*(u**2.0d0+v**2.0d0) )
    a          = sqrt(gama*p/rho)
    q_vel = sqrt(u**2.0d0 + v**2.0d0)
    if( (q_vel/a) <= 1.0d0 ) then
        Q_barra(i,j,1) = metric_jacobian(i,j)*(Q_barra(i-1,j,1)/metric_jacobian(i-1,j))
        Q_barra(i,j,2) = metric_jacobian(i,j)*(Q_barra(i-1,j,2)/metric_jacobian(i-1,j))
        Q_barra(i,j,3) = metric_jacobian(i,j)*(Q_barra(i-1,j,3)/metric_jacobian(i-1,j))
        p              = p_total/3.0d0
        u   = Q_barra(i-1,j,2)/Q_barra(i-1,j,1)
        v   = Q_barra(i-1,j,3)/Q_barra(i-1,j,1)
        Q_barra(i,j,4) = metric_jacobian(i,j)*(p/(gama-1.0d0) &
                        + (Q_barra(i,j,1)/(2.0d0*metric_jacobian(i,j)))*(u**2.0d0 + v**2.0d0))
    else
        Q_barra(i,j,1) = metric_jacobian(i,j)*(Q_barra(i-1,j,1)/metric_jacobian(i-1,j))
        Q_barra(i,j,2) = metric_jacobian(i,j)*(Q_barra(i-1,j,2)/metric_jacobian(i-1,j))
        Q_barra(i,j,3) = metric_jacobian(i,j)*(Q_barra(i-1,j,3)/metric_jacobian(i-1,j))
        Q_barra(i,j,4) = metric_jacobian(i,j)*(Q_barra(i-1,j,4)/metric_jacobian(i-1,j))
    end if
end do

! symmetry boundary

j = jmax
do i = 1, imax
    Q_barra(i,j,1) =  metric_jacobian(i,j)*(Q_barra(i,j-2,1)/metric_jacobian(i,j-2))
    Q_barra(i,j,2) =  metric_jacobian(i,j)*(Q_barra(i,j-2,2)/metric_jacobian(i,j-2))
    Q_barra(i,j,3) = -metric_jacobian(i,j)*(Q_barra(i,j-2,3)/metric_jacobian(i,j-2))
    Q_barra(i,j,4) =  metric_jacobian(i,j)*(Q_barra(i,j-2,4)/metric_jacobian(i,j-2))
end do

! do j = 1, jmax
!     do i = 1, imax
!         if (p(i,j) < 0.00003d0) then 
!                 print *, i,j,iter,p(i,j)
!         end if
!     end do     
! end do

end subroutine boundary_conditions_curv
