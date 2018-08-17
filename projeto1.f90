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
    dim = 4
    pi = dacos(-1.0d0)
    theta = 0.0d0*(180.0d0/pi)
    gama = 1.4d0
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
! add one more point on the index j due to the symmetry line
!
jmax = jmax + 1
call allocate_vars
!
!
!
call mesh
call metric_terms
call initial_conditions
!call output
call boundary_conditions
call output
! coordinate transformation
!
! loop - update solution
!
!
call deallocate_vars
!
close(1)
close(2)
!
end program proj1
