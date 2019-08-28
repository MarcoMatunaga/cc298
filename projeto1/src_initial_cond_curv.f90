subroutine initial_conditions_curv
    use vars
    implicit none
! change: but the definitions of boundary conditions in a generic way
real(8), dimension(:,:,:), allocatable          :: Q
real(8)                                         :: u, v, p, T, rho, e, a
real(8)                                         :: u_inf, v_inf, p_inf, rho_inf, a_inf
real(8)                                         :: u_til, v_til, p_til, rho_til, e_til


! inlet boundary

! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K
! free stream values
T = T_total
p_inf = p_total
rho_inf = p_inf/(R*T)
u_inf = 0.0d0
v_inf = 0.0d0
Q_barra         = 0.0d0
U_contravariant = 0.0d0
V_contravariant = 0.0d0
e = p_inf/(gama-1.0d0) + (rho_inf/2.0d0)*(u_inf**2.0d0 + v_inf**2.0d0)
! speed of sound at the freestream
a_inf = sqrt(gama*p_inf/rho_inf)

! calculate non-dimensional value
rho_til = 1.0d0
e_til = e/(rho_inf*a_inf**2.0d0)
p_til = p_inf/(rho_inf*a_inf**2.0d0)
u_til = u_inf/a_inf
v_til = v_inf/a_inf

! Euler vectors in cartesian coordinates
! non-dimensional variables
allocate(Q(imax,jmax,dim))

do j = 1, jmax
    do i = 1, imax
        Q(i,j,1) = rho_til    
        Q(i,j,2) = rho_til*u_til
        Q(i,j,3) = rho_til*v_til
        Q(i,j,4) = e_til
    end do
end do

! coordinates transformation
do j = 1, jmax
    do i = 1, imax
        Q_barra(i,j,1) = metric_jacobian(i,j)*Q(i,j,1)
        Q_barra(i,j,2) = metric_jacobian(i,j)*Q(i,j,2)
        Q_barra(i,j,3) = metric_jacobian(i,j)*Q(i,j,3)
        Q_barra(i,j,4) = metric_jacobian(i,j)*Q(i,j,4)
    end do
end do  

! change the outlet values
p = p_total/3.0d0
rho = p/(R*T)
u = 0.0d0
v = 0.0d0
U_contravariant = 0.0d0
V_contravariant = 0.0d0
e = p/(gama-1.0d0) + (rho/2.0d0)*(u**2.0d0 + v**2.0d0)
! speed of sound at the outlet boundary
a = sqrt(gama*p/rho)

! calculate non-dimensional value
rho_til = rho/rho_inf
e_til = e/(rho_inf*a_inf**2.0d0)
p_til = p/(rho_inf*a_inf**2.0d0)
u_til = u/a_inf
v_til = v/a_inf

i = imax
do j = 1, jmax
    Q(i,j,1)       = rho_til
    Q(i,j,4)       = e_til
    Q_barra(i,j,1) = metric_jacobian(i,j)*Q(i,j,1)
    Q_barra(i,j,4) = metric_jacobian(i,j)*Q(i,j,4)
end do

do j = 1, jmax
    do i = 1, imax
        U_contravariant(i,j) = u_til*ksi_x(i,j) + v_til*ksi_y(i,j)
        V_contravariant(i,j) = u_til*eta_x(i,j) + v_til*eta_y(i,j)
    end do 
end do

deallocate(Q)

end subroutine initial_conditions_curv
