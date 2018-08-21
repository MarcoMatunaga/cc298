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
    CFL = 0.001d0
    !
    ! inicializar com um deltat - duvida
    !
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
    iter = 0
    max_iter = 100000
!
! add one more point on the index j due to the symmetry line
!
jmax = jmax + 1
call allocate_vars
!
call mesh
call metric_terms
call initial_conditions_curv
!call output
call boundary_conditions_curv
max_residue = 1.0d0
!
!
do while ( max_residue > 10e-9 .and. iter < max_iter)
    call fluxes_curvilinear
    do j = 1, jmax
            do i = 1, imax
            a(i,j) = sqrt(gama*p(i,j)*metric_jacobian(i,j)/Q_barra(i,j,1))
            delta_t(i,j) = CFL/(max( abs(U_contravariant(i,j)) + a(i,j)*sqrt(ksi_x(i,j)**2.0d0 + ksi_y(i,j)**2.0d0 ), &
                             abs(V_contravariant(i,j)) + a(i,j)*sqrt(eta_x(i,j)**2.0d0 + eta_y(i,j)**2.0d0 )))
            end do
    end do
    call euler_explicit
    call boundary_conditions_curv
    iter = iter + 1
end do
call output
!
!
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
