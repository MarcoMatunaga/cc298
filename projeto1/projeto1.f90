!
! program project 1 cc298
!
program proj1
        use vars 
        use output_routines
        implicit none
integer(4)              :: i_cfl
integer(4)                                   :: ramp, CFL_ramp_size
real(8),dimension(:,:),allocatable           :: CFL_ramp

allocate(CFL_ramp(CFL_ramp_size,2))
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
    !
    ! choose the time marching method method
    ! time_method = 1 euler explicit 
    ! time_method = 2 euler implicit (beam warming)
    !
    time_method = 2
    !
    ! choose dissipation (explicit)
    ! which_diss = 1 D4
    ! which_diss = 2 D2
    ! which_diss = 3 nonlinear
    !
    which_diss = 1
    !eps_e = 10.00d0 ! good value for the artificial dissipation of second differences
    !eps_e = 5.0d0  ! good value for the artificial dissipation of fourth differences
    !
    ! set dissipation parameters
    !
    eps_dis_e = 5.0d0
    eps_dis_i = 4.0d0*eps_dis_e
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
    nsave = 0
    iter = 0
    max_iter = 25000
    a_cr = sqrt((2.0d0*gama)*((gama-1.0d0)/(gama+1.0d0))*c_v*T_total)
    !
    ! max CFL if the CFL_ramp is on
    ! ramp = 1 there is a CFL ramp
    ! ramp = 0 no
    ramp = 1
    !
    !
    CFL_ramp_size = 3
    CFL = 20.00d0
    ! ***********
    ! outra possibilidade subir o CFL 
    ! a cada determinado numero delta de iteracoes
    ! ***********
    !
    ! CFL_ramp(i,j) = iter, CFL
    !
    CFL_ramp(1,1) = 0.0d0
    CFL_ramp(2,1) = 100.0d0
    ! deal this as percentage
    CFL_ramp(3,1) = DBLE((3/100)*max_iter)
    !
    ! set the CFL numbers
    !
    CFL_ramp(1,2) = 0.150d0
    CFL_ramp(2,2) = 15.0d0
    CFL_ramp(3,2) = CFL
!
! add one more point on the index j due to the symmetry line
!
jmax = jmax + 1
call allocate_vars
!
call mesh
call metric_terms
!
!
call output_metric_terms
!
!
call initial_conditions_curv
call boundary_conditions_curv
call output_inicial
max_residue = 1.0d0
i_cfl = 1
!
!
do while ( max_residue > -15.0d0 .and. iter < max_iter )
    call fluxes_curvilinear 
    call output_fluxes
    do j = 1, jmax
            do i = 1, imax
               a(i,j) = sqrt(gama*p(i,j)*metric_jacobian(i,j)/Q_barra(i,j,1))
                if (ramp == 1 .and. iter == int(CFL_ramp(i_cfl,1)) ) then                
                    CFL = CFL_ramp(i_cfl,2)
                    i_cfl = i_cfl + 1
                end if
                delta_t(i,j) = CFL/(max( abs(U_contravariant(i,j)) + a(i,j)*sqrt(ksi_x(i,j)**2.0d0 + ksi_y(i,j)**2.0d0 ), &
                                       abs(V_contravariant(i,j)) + a(i,j)*sqrt(eta_x(i,j)**2.0d0 + eta_y(i,j)**2.0d0 )))
            end do
    end do
    !
    ! time marching
    !
    if (time_method == 1) call euler_explicit
    if (time_method == 2) call implicit_beam_warming
    !
    !
    if ( mod(iter,(max_iter/10)) == 0 ) then
         nsave = nsave + 1
         call output_tecplot
    end if
    iter = iter + 1
    write(*,*) iter, CFL
    !
    !
    call output_residue
    call boundary_conditions_curv
end do
!
! close the archive used to write the residue
!
close(5)
!
!
call output_final
call deallocate_vars
deallocate(CFL_ramp)
!
!
close(1)
close(2)
!
!
end program proj1
