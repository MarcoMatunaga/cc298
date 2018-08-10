!
! program project 1 cc298
!
program proj1
use vars 
implicit none
!
!
!
gama = 1.4d0
pi = 3.1415d0
delta_eta = 1.0d0
delta_ksi = 1.0d0
!
! here we consider the value of cv as (5/2)*R
! R is the universal perfect gas constant
R = 287.053d0
c_v = 2.5d0*R
!
!
!
open(1,file='mesh') 
read(1,*) imax,jmax,kmax
!
! add one more point on the index j due to the symmetry line
!
jmax = jmax + 1
allocate(meshx(imax,jmax), meshy(imax,jmax), Q(imax,jmax,4))
allocate(p(imax,jmax), T(imax,jmax), rho(imax,jmax) )
allocate(u(imax,jmax), v(imax,jmax), a(imax,jmax))
allocate(y_ksi(imax,jmax), y_eta(imax,jmax))
allocate(x_ksi(imax,jmax), x_eta(imax,jmax))
!
!
!
!T_total = 0.555556d0*531.2d0
T_total = 294.8d0
!p_total = 47.880258888889d0*2117.0d0
p_total = 101360.0d0
call mesh
call metric_terms
call initial_conditions
call boundary_conditions
!
!
deallocate(meshx, meshy, Q)
deallocate(p, T, rho)
deallocate(u, v, a)
deallocate(y_ksi, y_eta)
deallocate(x_ksi, x_eta)
close(1)
close(2)
!
!
!
end program proj1
