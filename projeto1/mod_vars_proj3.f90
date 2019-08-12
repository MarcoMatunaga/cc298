module vars_proj3
    use vars
implicit none
    real(8)                      :: Height, Length
    real(8)                      :: shock_angle, mach_inlet, turn_angle, mach_inlet_n
    real(8)                      :: mach_2, mach_2_n
    real(8)                      :: x_wall, x_a, y_sol, shock_r_angle, lambda_shock, khi, x_3, z, x_ana
contains
    
    subroutine analytical
        implicit none
        real(8)                         :: p, rho

        shock_angle  = shock_angle*(3.1415d0/180.0d0)
        
        ! calculo para determinar as condicoes pos choque
        
        turn_angle   = datan( (2.0d0/tan(shock_angle)) * ( (mach_inlet*sin(shock_angle))**2.0d0 &
                       - 1.0d0 )/( (mach_inlet**2.0d0)*(gama + cos(2.0d0*shock_angle)) + 2.0d0)  )
        mach_inlet_n = mach_inlet*sin(shock_angle) 
        mach_2_n     = sqrt((2.0d0 + (gama - 1.0d0)*mach_inlet_n**2.0d0)/( 2.0d0*gama*mach_inlet_n**2.0d0 - (gama - 1.0d0) ) )
        mach_2       =  mach_2_n/sin(shock_angle - turn_angle)
        
        ! Determination of the point of reflection 
        ! in the bottom wall

        x_wall = Height/tan(shock_angle)
        x_a = tan(3.1415d0/2.0d0 - shock_angle)*(Height - y_sol)

        lambda_shock  = sqrt( (mach_2**2.0d0 - 1.0d0)**2.0d0 - 3.0d0*(1.0d0 + (gama-1.0d0)*mach_2**2.0d0/2.0d0 )*(1.0d0 &
                        + (gama+1.0d0)*mach_2**2.0d0/2.0d0)*tan(turn_angle)**2.0d0 ) 
        khi           = ( (mach_2**2.0d0 - 1.0d0)**3.0d0 - 9.0d0*(1.0d0 + (gama-1.0d0)*mach_2**2.0d0/2.0d0 )*(1.0d0 &
                        + (gama-1.0d0)*mach_2**2.0d0/2.0d0 &
                        + (gama+1.0d0)*mach_2**4.0d0/4.0d0 )*tan(turn_angle)**2.0d0)/(lambda_shock**3.0d0)
        shock_r_angle = atan( (mach_2**2.0d0 - 1.0d0 + 2.0d0*lambda_shock*cos((4.0d0*3.1415d0 + acos(khi))/3.0d0) )/(3.0d0*(1.0d0&
                        +(gama-1.0d0)*(mach_2**2.0d0)/2.0d0)*tan(turn_angle)) )
        write(*,*) shock_r_angle*(180.0d0/dacos(-1.0d0)), mach_2
        z = y_sol/sin(shock_r_angle - turn_angle)
        x_3 = z*cos(shock_r_angle - turn_angle) + x_wall

        ! Create the output solution file

        open(9,file = 'analytical.dat')

        ! Values for region 1 or initial conditions

        x_ana = 0.0d0
            do while (x_ana < x_a)
                 p  = 1.0d0/gama
                 rho = 1.0d0
                 write(9,'(3es11.3e2)') x_ana, p, rho
                 x_ana = x_ana + 0.01d0  
            end do
            
            ! Values for region 2 or post-shock solution conditions
            
            do while (x_ana >= x_a .and. x_ana < x_3)
                 rho =  (gama + 1.0d0)*mach_inlet_n**2.0d0/(2.0d0 + (gama - 1.0d0)*mach_inlet_n**2.0d0)
                 p   =  (1.0d0 + (2.0d0*gama/( gama + 1.0d0 ))*(mach_inlet_n**2.0d0 - 1.0d0) ) / gama
                 write(9,'(3es11.3e2)') x_ana, p, rho
                 x_ana = x_ana + 0.01d0  
            end do
            
            ! Values for region 3 or post-reflected shock conditions
            
            do while (x_ana >= x_3 .and. x_ana < Length)
                 mach_2_n = mach_2*sin(shock_r_angle)
                 rho      = (gama + 1.0d0)*mach_inlet_n**2.0d0/(2.0d0 + (gama - 1.0d0)*mach_inlet_n**2.0d0) * (gama &
                            + 1.0d0)*mach_2_n**2.0d0/(2.0d0 + (gama - 1.0d0)*mach_2_n**2.0d0)
                 p     =  ((1.0d0 + (2.0d0*gama/( gama + 1.0d0 ))*(mach_inlet_n**2.0d0 - 1.0d0) ) / gama) * ((1.0d0 &
                            + (2.0d0*gama/( gama + 1.0d0 ))*(mach_2_n**2.0d0 - 1.0d0) ) )

                 write(9,'(3es11.3e2)') x_ana, p, rho

                 x_ana = x_ana + 0.01d0  
            end do

        close(9)

    end subroutine analytical

end module vars_proj3