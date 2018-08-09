!
! implementation of boundary conditions
!
subroutine boundary_conditions
use vars
implicit none
! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K
!
! inlet boundary
! 
theta = 0.0d0*(180.0d0/pi)
i = 1
do j = 1, jmax - 1
	u(i,j) = Q(i,j,2)/Q(i,j,1) 
	v(i,j) = u(i,j)*tan(theta)		
	a(i,j) = sqrt(2.0d0*gama*((gama-1.0d0)/(gama+1.0d0))*c_v*T_total)
	T(i,j) = T_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*(1.0d0+tan(theta)**2.0)*((u(i,j)/a(i,j))**2.0d0))
	p(i,j) = p_total*(1.0d0-((gama-1.0d0)/(gama+1.0d0))*(1.0d0+tan(theta)**2.0)*((u(i,j)/a(i,j))**2.0d0))**(gama/(gama-1.0d0))
	Q(i,j,4) = p(i,j)/(gama-1.0d0) + (rho(i,j)/2.0d0)*(u(i,j)**2 + v(i,j)**2)
end do
!
! outlet boundary
!
i = imax
do j = 1, jmax - 1 
		p(i,j) = p_total/3.0d0
                Q(i,j,4) = p(i,j)/(gama-1.0d0) + (rho(i,j)/2.0d0)*(u(i,j)**2 + v(i,j)**2)
end do
!
!
!
q_vel = u**2 +v**2
!
!
!
end subroutine boundary_conditions
