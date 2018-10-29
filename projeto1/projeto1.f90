!
! program project 1 cc298
!
program proj1
        use vars 
        use output_routines
        implicit none
integer(4)                                   :: i_cfl
integer(4)                                   :: ramp, CFL_ramp_size
real(8),dimension(:),allocatable             :: CFL_ramp
integer(4),dimension(:),allocatable          :: CFL_ramp_iter
!
!
namelist /PAR_Flow/ gama, R, T_total,p_total 
namelist /PAR_Time/ max_iter,time_method,res_conv
namelist /PAR_Numeric/ which_diss, k2, k4, eps_dis_e, dis_factor, ramp, CFL_ramp_size, CFL
namelist /PAR_Others/ dim 
!
!
open(8,file='inputs')
read(8,PAR_Flow)
read(8,PAR_Time)
read(8,PAR_Numeric)
read(8,PAR_Others)
close(8)
!
!
allocate(CFL_ramp(CFL_ramp_size),CFL_ramp_iter(CFL_ramp_size))
!
!
open(1,file='mesh') 
read(1,*) imax,jmax
!
!
!
    pi = dacos(-1.0d0)
    theta = 0.0d0*(180.0d0/pi)
    delta_eta = 1.0d0
    delta_ksi = 1.0d0
    ! 
    !
    ! choose the time marching method method
    ! time_method = 1 euler explicit 
    ! time_method = 2 euler implicit (beam warming)
    !
    !
    ! choose dissipation (explicit)
    ! which_diss = 1 D4
    ! which_diss = 2 D2
    ! which_diss = 3 nonlinear
    !
    !eps_e = 10.00d0 ! good value for the artificial dissipation of second differences
    !eps_e = 5.0d0  ! good value for the artificial dissipation of fourth differences
    !
    ! set dissipation parameters
    !
    eps_dis_i = dis_factor*eps_dis_e
    !
    ! here we consider the value of cv as (5/2)*R
    ! R is the universal perfect gas constant
    !
    c_v = 2.5d0*R
    !T_total = 0.555556d0*531.2d0
    !p_total = 47.880258888889d0*2117.0d0
    nsave = 0
    iter = 0
    a_cr = sqrt((2.0d0*gama)*((gama-1.0d0)/(gama+1.0d0))*c_v*T_total)
    !
    ! max CFL if the CFL_ramp is on
    ! ramp = 1 there is a CFL ramp
    ! ramp = 0 no
    !
    !
    ! ***********
    ! outra possibilidade subir o CFL 
    ! a cada determinado numero delta de iteracoes
    ! ***********
    !
    ! CFL_ramp = iter, CFL
    !
    CFL_ramp_iter(1) = 0
    CFL_ramp_iter(2) = 200
    ! deal this as percentage
    CFL_ramp_iter(3) = 750
    !
    ! set the CFL numbers
    !
    CFL_ramp(1) = 5.00d0
    CFL_ramp(2) = 10.0d0
    CFL_ramp(3) = CFL
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
do while ( max_residue > res_conv .and. iter < max_iter )
    call fluxes_curvilinear 
    call output_fluxes
    do j = 1, jmax
            do i = 1, imax

                a(i,j) = sqrt(gama*p(i,j)*metric_jacobian(i,j)/Q_barra(i,j,1))
                if (ramp == 1 .and. iter == CFL_ramp_iter(i_cfl) ) then                
                    CFL = CFL_ramp(i_cfl)
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
    if (time_method == 3) call pulliam_chausse
    ! if (time_method == 4) call sw_1st
    ! if (time_method == 5) call sw_2nd
    !
    !
    if ( mod(iter,(max_iter/10)) == 0 ) then
         nsave = nsave + 1
         call output_tecplot
    end if
    iter = iter + 1
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
deallocate(CFL_ramp,CFL_ramp_iter)
!
!
close(1)
close(2)
!
!
end program proj1
