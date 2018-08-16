!
! program project 1 cc298
!
program proj1
use vars 
implicit none
!
!
!
theta = 0.0d0*(180.0d0/pi)
gama = 1.4d0
pi = 3.1415d0
delta_eta = 1.0d0
delta_ksi = 1.0d0
!
! here we consider the value of cv as (5/2)*R
! R is the universal perfect gas constant
!
R = 287.053d0
c_v = 2.5d0*R
!T_total = 0.555556d0*531.2d0
T_total = 294.8d0
!p_total = 47.880258888889d0*2117.0d0
p_total = 101360.0d0
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
allocate(Q(imax,jmax,4))
allocate(phi_jacobian(imax,jmax), theta_jacobian(imax,jmax))
allocate(a1_jacobian(imax,jmax), metric_jacobian(imax,jmax))
allocate(p(imax,jmax), T(imax,jmax))
allocate(u(imax,jmax), v(imax,jmax), a(imax,jmax), q_vel(imax,jmax))
allocate(a_cr(imax,jmax))
allocate(y_ksi(imax,jmax), y_eta(imax,jmax))
allocate(x_ksi(imax,jmax), x_eta(imax,jmax))
allocate(ksi_x(imax,jmax), ksi_y(imax,jmax))
allocate(eta_x(imax,jmax), eta_y(imax,jmax))
allocate(A_barra(imax,jmax), B_barra(imax,jmax), M_barra(imax,jmax))
allocate(E_barra(imax,jmax,4), F_barra(imax,jmax,4))
!
!
!
call mesh
call metric_terms
call output
call initial_conditions
call boundary_conditions
!
!
deallocate(meshx, meshy, Q)
deallocate(p, T)
deallocate(u, v, a, q_vel, a_cr)
deallocate(y_ksi, y_eta)
deallocate(x_ksi, x_eta)
deallocate(ksi_x, ksi_y)
deallocate(eta_x, eta_y)
deallocate(phi_jacobian, theta_jacobian, a1_jacobian)
deallocate(metric_jacobian)
deallocate(A_barra, B_barra, M_barra)
deallocate(E_barra, F_barra)
!
!
close(1)
close(2)
!
end program proj1
