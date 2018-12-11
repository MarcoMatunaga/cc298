
! program project 1 cc298

program proj1
        use vars 
        use vars_proj3
        use vars_sw
        use output_routines
        use fluxes_pos_neg
        implicit none
        real(8)                         :: start, finish
        real(8)                         :: start_iter, finish_iter
        real(8)                         :: start_save, finish_save, t_save
!****
! create a fortran program with command line inputs
!****
namelist /PAR_Flow/ gama, R, T_total, p_total 
namelist /PAR_Flow_shock/ mach_inlet, shock_angle
namelist /PAR_Time/ max_iter,time_method,res_conv
namelist /PAR_Numeric/ which_diss, k2, k4, eps_dis_e, dis_factor, ramp, CFL_ramp_size, CFL
namelist /PAR_Analytical/ y_sol
namelist /PAR_Others/ dim, which_boundary, imax, jmax, Height, Length

open(8,file='inputs')
read(8,PAR_Flow)
read(8,PAR_Flow_shock)
read(8,PAR_Time)
read(8,PAR_Numeric)
read(8,PAR_Analytical)
read(8,PAR_Others)
close(8)

call allocate_vars

    delta_eta = 1.0d0
    delta_ksi = 1.0d0

if (which_boundary == 1) then 
    open(1,file='mesh') 
    read(1,*) imax,jmax
    jmax = jmax + 1
    call mesh
    call metric_terms
    call output_metric_terms
else 
    call mesh_proj3
    call metric_terms
    call output_metric_terms
    call analytical
end if

    pi = dacos(-1.0d0)
    theta = 0.0d0*(180.0d0/pi)
    ! choose the time marching method method
    ! time_method = 1 euler explicit 
    ! time_method = 2 euler implicit (beam warming)
    
    ! choose dissipation (explicit)
    ! which_diss = 1 D4
    ! which_diss = 2 D2
    ! which_diss = 3 nonlinear
    
    !eps_e = 10.00d0 ! good value for the artificial dissipation of second differences
    !eps_e = 5.0d0  ! good value for the artificial dissipation of fourth differences
    
    ! set dissipation parameters
    
    eps_dis_i = dis_factor*eps_dis_e
    
    ! here we consider the value of cv as (5/2)*R
    ! R is the universal perfect gas constant
    
    c_v = 2.5d0*R

    !T_total = 0.555556d0*531.2d0
    !p_total = 47.880258888889d0*2117.0d0
    nsave = 0
    t_save = 0.0d0
    iter = 0
    a_cr = sqrt((2.0d0*gama)*((gama-1.0d0)/(gama+1.0d0))*c_v*T_total)
    
! add one more point on the index j due to the symmetry line

if (which_boundary == 1) then 
    call initial_conditions_curv
    call boundary_conditions_curv
else
    call initial_cond_proj3
    call boundary_proj3
end if 

max_residue = 1.0d0
call output_inicial 

    call cpu_time(start)

do 
    if ( max_residue < res_conv .or. iter > max_iter ) exit

    if (time_method == 1 .or. time_method == 2 .or. time_method == 3) call fluxes_curvilinear 
    if (time_method == 4 .or. time_method == 5 .or. time_method == 11 .or. time_method == 12) then
        call allocate_vars_sw
        call calculate_fluxes(E_pos,E_neg,F_pos,F_neg)
    end if 
    call output_fluxes
    call calculate_CFL
    
    ! time marching
        call cpu_time(start_iter)
    if (time_method == 1) call euler_explicit
    if (time_method == 2) call implicit_beam_warming
    if (time_method == 3) call pulliam_chausse_block
    if (time_method == 4) call sw_1st
    if (time_method == 5) call sw_2nd
    ! if (time_method == 6) call ausm_plus_1st
    ! if (time_method == 7) call ausm_plus_2nd
    ! if (time_method == 8) call van_leer_1st
    ! if (time_method == 9) call van_leer_2nd
    if (time_method == 10) call pulliam_chausse
    if (time_method == 11) call sw_1st_a
    if (time_method == 12) call sw_2nd_a
        call cpu_time(finish_iter)
        print '(i0,"Time per iterration = ",f6.3," seconds.")',iter,finish_iter-start_iter

    iter = iter + 1
    if ( mod(iter,(max_iter/10)) == 0 ) then
         call cpu_time(start_save)
         nsave = nsave + 1
         call output_tecplot
         call cpu_time(finish_save)
         t_save = t_save + finish_save - start_save
    end if
    
    call output_residue
    if (which_boundary == 1) call boundary_conditions_curv
    if (which_boundary == 2) call boundary_proj3

end do

        call cpu_time(finish)
        print '("Time = ",f6.3," seconds.")',finish-start - t_save
! close the archive used to write the residue

close(5)

call output_final
call deallocate_vars

close(2)

end program proj1
