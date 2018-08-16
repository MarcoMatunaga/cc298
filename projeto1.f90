!
! program project 1 cc298
!
program proj1
use vars 
implicit none
!
!
!
open(1,file='mesh') 
read(1,*) imax,jmax
!
! add one more point on the index j due to the symmetry line
!
jmax = jmax + 1
allocate(meshx(imax,jmax), meshy(imax,jmax))
allocate(Q(imax,jmax,dim))
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
allocate(E_barra(imax,jmax,dim), F_barra(imax,jmax,dim))
!
!
!
call mesh
call metric_terms
call output
call initial_conditions
call boundary_conditions
! coordinate transformation
!
! loop - update solution
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
