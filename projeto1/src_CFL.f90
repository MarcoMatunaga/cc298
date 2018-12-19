subroutine calculate_CFL
    use vars
    implicit none
    real(8)                                      :: a, p, rho
    integer(4)                                   :: i_cfl
    real(8),dimension(:),allocatable             :: CFL_ramp
    integer(4),dimension(:),allocatable          :: CFL_ramp_iter
    
    allocate(CFL_ramp(CFL_ramp_size),CFL_ramp_iter(CFL_ramp_size))

    ! max CFL if the CFL_ramp is on
    ! ramp = 1 there is a CFL ramp
    ! ramp = 0 no CFL ramp
    
    
    ! ***********
    ! outra possibilidade subir o CFL 
    ! a cada determinado numero delta de iteracoes
    ! ***********

    i_cfl = 1

    ! CFL_ramp = iter, CFL
    
    CFL_ramp_iter(1) = 0
    CFL_ramp_iter(2) = 200

    ! deal this as percentage
    CFL_ramp_iter(3) = 750
    
    ! set the CFL numbers
    
    CFL_ramp(1) = 5.00d0
    CFL_ramp(2) = 10.0d0
    CFL_ramp(3) = CFL

    do j = 1, jmax
        do i = 1, imax
            rho = Q_barra(i,j,1)/metric_jacobian(i,j)
            p = (gama-1.0d0) * (Q_barra(i,j,4)/metric_jacobian(i,j) & 
               - 0.5d0*( (Q_barra(i,j,2)/metric_jacobian(i,j))**2.0d0 &
               + (Q_barra(i,j,3)/metric_jacobian(i,j))**2.0d0)/rho)
            a  = sqrt(gama*p*metric_jacobian(i,j)/Q_barra(i,j,1))
            if (ramp == 1 .and. iter == CFL_ramp_iter(i_cfl) ) then                
                    CFL = CFL_ramp(i_cfl)
                    i_cfl = i_cfl + 1
            end if
            delta_t(i,j) = CFL/(max( abs(U_contravariant(i,j)) + a*sqrt(ksi_x(i,j)**2.0d0 + ksi_y(i,j)**2.0d0 ), &
                                     abs(V_contravariant(i,j)) + a*sqrt(eta_x(i,j)**2.0d0 + eta_y(i,j)**2.0d0 )))
        end do
    end do

    deallocate(CFL_ramp,CFL_ramp_iter)

end subroutine calculate_CFL