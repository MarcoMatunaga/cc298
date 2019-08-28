
! implementation of boundary conditions

subroutine boundary_conditions_curv
    use vars
    implicit none
        real(8)                         :: u, v, p, rho, T, a, q_vel, en
        real(8)                         :: p_inf, rho_inf, a_inf, e_inf

! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K

! freestream variables
p_inf = p_total
rho_inf = p_inf/(R*T_total)
a_inf = sqrt(gama*p_inf/rho_inf) 

! inlet boundary
i = 1
do j = 2, jmax - 1

    u  = (Q_barra(i+1,j,2)/Q_barra(i+1,j,1))/a_inf
    v  = u*dtan(theta)/a_inf  
    T  = T_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*((u/a_cr)**2.0d0))
    p  = (p_total/(rho_inf*a_inf**2.0d0))*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*((u/a_cr)**2.0d0))**(gama/(gama-1.0d0))
    rho = (p/(R*T))/rho_inf
    Q_barra(i,j,1) = metric_jacobian(i,j)*rho
    Q_barra(i,j,2) = metric_jacobian(i,j)*rho*u
    Q_barra(i,j,3) = metric_jacobian(i,j)*rho*v
    Q_barra(i,j,4) = metric_jacobian(i,j)*(p/(gama-1.0d0) &
                     + (0.50d0*rho)*(u**2.0d0 + v**2.0d0))
    
    U_contravariant(i,j) = u*ksi_x(i,j) + v*ksi_y(i,j)
    V_contravariant(i,j) = u*eta_x(i,j) + v*eta_y(i,j)

end do

! outlet boundary

i = imax
do j = 2, jmax - 1

    ! Mach calculation
    rho = Q_barra(i,j,1)/metric_jacobian(i,j)
    u   = Q_barra(i,j,2)/Q_barra(i,j,1)
    v   = Q_barra(i,j,3)/Q_barra(i,j,1)
    p = (gama-1.0d0) * (Q_barra(i,j,4)/metric_jacobian(i,j) & 
               - 0.5d0*rho*(u**2.0d0+v**2.0d0) )
    a          = sqrt(gama*p/rho)
    q_vel = sqrt(u**2.0d0 + v**2.0d0)

    ! extrapolate the properties
    rho = Q_barra(i-1,j,1)/metric_jacobian(i-1,j)
    u   = Q_barra(i-1,j,2)/Q_barra(i-1,j,1)
    v   = Q_barra(i-1,j,3)/Q_barra(i-1,j,1)

    if ( (q_vel/a) <= 1.0d0) then
        p   = (p_total/3.0d0)/(rho_inf*a_inf**2.0d0)
    else
        en  = Q_barra(i-1,j,4)/metric_jacobian(i-1,j)
        p   = (gama - 1.0d0)*(en - 0.50d0*rho*(u**2.0d0+v**2.0d0))*(rho_inf*a_inf**2.0d0)
    endif

    Q_barra(i,j,1) = metric_jacobian(i,j)*rho
    Q_barra(i,j,2) = metric_jacobian(i,j)*rho*u
    Q_barra(i,j,3) = metric_jacobian(i,j)*rho*v
    Q_barra(i,j,4) = metric_jacobian(i,j)*(p/(gama-1.0d0) &
                     + (0.50d0*rho*(u**2.0d0 + v**2.0d0)))

    U_contravariant(i,j) = u*ksi_x(i,j) + v*ksi_y(i,j)
    V_contravariant(i,j) = u*eta_x(i,j) + v*eta_y(i,j)

end do

! symmetry boundary

j = jmax
do i = 1, imax
   
    rho = Q_barra(i,j-2,1)/metric_jacobian(i,j-2)
    u   = Q_barra(i,j-2,2)/Q_barra(i,j-2,1)
    v   = -Q_barra(i,j-2,3)/Q_barra(i,j-2,1)
    en  = Q_barra(i,j-2,4)/metric_jacobian(i,j-2)
    p   = (gama - 1.0d0)*(en - 0.50d0*rho*(u**2.0d0+v**2.0d0))
    
    Q_barra(i,j,1) = metric_jacobian(i,j)*rho
    Q_barra(i,j,2) = metric_jacobian(i,j)*rho*u
    Q_barra(i,j,3) = metric_jacobian(i,j)*rho*v
    Q_barra(i,j,4) = metric_jacobian(i,j)*(p/(gama-1.0d0) &
                     + (0.50d0*rho*(u**2.0d0 + v**2.0d0)))

    U_contravariant(i,j) = u*ksi_x(i,j) + v*ksi_y(i,j)
    V_contravariant(i,j) = u*eta_x(i,j) + v*eta_y(i,j)

end do

! lower boundary (wall boundary)

j = 1
do i = 1, imax
       
        ! u e v sao determinados pela matriz jacobiana de transformacao
        !so lebrar dos termos contrvariante das variaveis

    rho = Q_barra(i,j+1,1)/metric_jacobian(i,j+1)
    u   = Q_barra(i,j+1,2)/Q_barra(i,j+1,1)
    v   = Q_barra(i,j+1,3)/Q_barra(i,j+1,1)
    en  = Q_barra(i,j+1,4)/metric_jacobian(i,j+1)
    p   = (gama - 1.0d0)*(en - 0.50d0*rho*(u**2.0d0+v**2.0d0)) 

    Q_barra(i,j,1) = metric_jacobian(i,j)*rho
    Q_barra(i,j,2) = metric_jacobian(i,j)*rho*u
    Q_barra(i,j,3) = metric_jacobian(i,j)*rho*v
    Q_barra(i,j,4) = metric_jacobian(i,j)*(p/(gama-1.0d0) &
                     + (0.50d0*rho*(u**2.0d0 + v**2.0d0)))

    U_contravariant(i,j) = u*ksi_x(i,j) + v*ksi_y(i,j)
    V_contravariant(i,j) = 0.0d0

end do

! do j = 1, jmax
!     do i = 1, imax
!         if (p(i,j) < 0.00003d0) then 
!                 print *, i,j,iter,p(i,j)
!         end if
!     end do     
! end do

end subroutine boundary_conditions_curv
