module vars
implicit none
! variaveis do fluido e do escoamento	
real(8),dimension(:,:),allocatable  		 :: T, p, u, v, rho, q_vel
real(8),dimension(:,:),allocatable 		 :: Q 
real(8)						 :: T_total, p_total
! constantes do fluido
real(8)						 :: gama, c_v
! constantes do escoamento 
real(8)						 :: theta 
! constantes matematicas
real(8)						 :: pi, dummy
! indices dos vetores --> geometria, malha
real(8),dimension(:,:),allocatable               :: meshx, meshy
integer(4)				         :: i,j,k
integer(4)					 :: imax,jmax,kmax
! variaveis numericas 
!
!
!
gama = 1.4d0
pi = 3.1415d0
c_v = 
!
!
!
end module vars
!
!
!
program proj1
use vars 
implicit none
!
!
!
open(1,file='mesh') 
read(1,*) imax,jmax,kmax
!
! add one more point on the index j due to the symmetry line
!
jmax = jmax + 1
allocate(meshx(imax,jmax), meshy(imax,jmax))
allocate(p(imax,jmax), T(imax,jmax), rho(imax,jmax) )
allocate(u(imax,jmax), v(imax,jmax))
!
!
T_total = 0.555556d0*531.2d0
p_total = 47.880258888889d0*2117.0d0
call mesh
call initial_conditions
call boundary_conditions
!
!
deallocate(meshx, meshy)
deallocate(p, T, rho)
deallocate(u, v)
!
!
!
end program proj1
!
!
!
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
rho =
!
! Euler vectors
!
do j = 1, jmax
	do i = 1, imax
		Q(i,j,1) = rho(i,j)	
		Q(i,j,2) = rho(i,j)*u(i,j)
		Q(i,j,3) = rho(i,j)*v(i,j)
		Q(i,j,4) = p(i,j)/(gama-1.0d0) + (rho/2.0d0)*(u(i,j)**2 + v(i,j)**2)
	end do
end do
!
!
!
end subroutine initial_conditions
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
	a(i,j) = sqrt(2.0d0*gama*( (gama - 1.0d0)/(gama + 1.0d0) )*c_v*T_total)
	T(i,j) = T_total*(1.0d0 - ( (gama-1.0d0)/(gama+1.0d0) )*(1.0d0+tan(theta)**2.0)*( u(i,j)/a(i,j) )**2.0d0) )
	p(i,j) = p_total*(1.0d0 - ( (gama-1.0d0)/(gama+1.0d0) )*(1.0d0+tan(theta)**2.0)*( u(i,j)/a(i,j) )**2.0d0) )**(gama/(gama-1.0d0))
	Q(i,j,4) = p(i,j)/(gama-1.0d0) + (rho/2.0d0)*(u(i,j)**2 + v(i,j)**2)
end do
!
! outlet boundary
!
i = imax
do j = 1, jmax - 1 
		p(i,j) = p_total/3.0d0
                Q(i,j,4) = p(i,j)/(gama-1.0d0) + (rho/2.0d0)*(u(i,j)**2 + v(i,j)**2)
end do
!
!
!
q_vel = u**2 +v**2
!
!
!
end subroutine boundary_conditions
!
!
!
subroutine implicit_beam_warming
use vars
implicit none

end subroutine implicit_beam_warming
!
!
!
subroutine mesh
use vars
implicit none
!
! construct the mesh 
!
do j = 1, jmax
	do i = 1, imax
		read(1,*) meshx(i,j),meshy(i,j), dummy ! dummy para ler a coordenada z
	end do
end do
close(1)
!
! arquivo tecplot
!
open(2,file='mesh.dat')
write(2,*) 'TITLE = "Projeto1" '
write(2,*) 'VARIABLES = "X" "Y" '
WRITE(2,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
do j = 1, jmax
	do i = 1, imax
		write(2,*) meshx(i,j), meshy(i,j)
	end do
end do
close(2)
!
!
!
end subroutine mesh
