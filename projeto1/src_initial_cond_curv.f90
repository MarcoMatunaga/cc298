subroutine initial_conditions_curv
    use vars
    implicit none

real(8), dimension(:,:,:), allocatable  :: Q
real(8)                                 :: u, v, p, T

! inlet boundary

! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K

T = T_total
p = p_total
u = 0.0d0
v = 0.0d0
Q_barra         = 0.0d0
U_contravariant = 0.0d0
V_contravariant = 0.0d0

! Euler vectors in cartesian coordinates

allocate(Q(imax,jmax,dim))

do j = 1, jmax
    do i = 1, imax
        Q(i,j,1) = p/(R*T)    
        Q(i,j,2) = Q(i,j,1)*u
        Q(i,j,3) = Q(i,j,1)*v
        Q(i,j,4) = p/(gama-1.0d0) + (Q(i,j,1)/2.0d0)*(u**2.0d0 + v**2.0d0)
    end do
end do

do j = 1, jmax
    do i = 1, imax
        Q_barra(i,j,1) = metric_jacobian(i,j)*Q(i,j,1)
        Q_barra(i,j,2) = metric_jacobian(i,j)*Q(i,j,2)
        Q_barra(i,j,3) = metric_jacobian(i,j)*Q(i,j,3)
        Q_barra(i,j,4) = metric_jacobian(i,j)*Q(i,j,4)
    end do
end do  

! change the outlet values

i = imax
do j = 1, jmax
    p              = p_total/3.0d0
    Q(i,j,1)       = p/(R*T)
    Q(i,j,4)       = p/(gama-1.0d0) + (Q(i,j,1)/2.0d0)*(u**2.0d0 + v**2.0d0)
    Q_barra(i,j,1) = metric_jacobian(i,j)*Q(i,j,1)
    Q_barra(i,j,4) = metric_jacobian(i,j)*Q(i,j,4)
end do
!
!
!
do j = 1, jmax
    do i = 1, imax
        U_contravariant(i,j) = u*ksi_x(i,j) + v*ksi_y(i,j)
        V_contravariant(i,j) = u*eta_x(i,j) + v*eta_y(i,j)
    end do 
end do
!
!
! 
deallocate(Q)
!
!
end subroutine initial_conditions_curv
