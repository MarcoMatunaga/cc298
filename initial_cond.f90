subroutine initial_conditions
use vars
implicit none
!
! inlet boundary
! 
! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K
!
T = T_total
p = p_total
u = 0.0d0
v = 0.0d0
!
! Euler vectors
!
do j = 1, jmax
    do i = 1, imax
        Q(i,j,1) = p(i,j)/(R*T(i,j))    
        Q(i,j,2) = Q(i,j,1)*u(i,j)
        Q(i,j,3) = Q(i,j,1)*v(i,j)
        Q(i,j,4) = p(i,j)/(gama-1.0d0) + (Q(i,j,1)/2.0d0)*(u(i,j)**2 + v(i,j)**2)
    end do
end do
!
! change the outlet values
!
i = imax
do j = 1, imax
    p(i,j) = p_total/3.0d0
end do
!
!
end subroutine initial_conditions
