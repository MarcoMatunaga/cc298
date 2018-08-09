subroutine initial_conditions
use vars
implicit none
!
! inlet boundary
! 
! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K
!
!
u = 0.0d0
v = 0.0d0
T = T_total
p = p_total
rho = p/(R*T)
!
! Euler vectors
!
do j = 1, jmax
	do i = 1, imax
		Q(i,j,1) = rho(i,j)	
		Q(i,j,2) = rho(i,j)*u(i,j)
		Q(i,j,3) = rho(i,j)*v(i,j)
		Q(i,j,4) = p(i,j)/(gama-1.0d0) + (rho(i,j)/2.0d0)*(u(i,j)**2 + v(i,j)**2)
	end do
end do
!
!
end subroutine initial_conditions
