!
! euler explicit
!
subroutine euler_explicit
        use vars
        implicit none
real(8)                                  :: eps_e
real(8),dimension(:,:,:), allocatable    :: Q_dis, D4_ksi, D4_eta
!
!
!
eps_e = 0.03d0
allocate(Q_dis(imax,jmax,dim), D4_ksi(imax,jmax,dim), D4_eta(imax,jmax,dim) )
!
! vector q_dis apenas para facilitar o calculo da dissipacao artificial
!
do j = 1, jmax
        do i = 1, imax
            Q_dis(i,j,1) = Q_barra(i,j,1)/metric_jacobian(i,j)
            Q_dis(i,j,2) = Q_barra(i,j,2)/metric_jacobian(i,j)
            Q_dis(i,j,3) = Q_barra(i,j,3)/metric_jacobian(i,j)
            Q_dis(i,j,4) = Q_barra(i,j,4)/metric_jacobian(i,j)
        end do
end do
!
!
!
do j = 3, jmax - 2
        do i = 3, imax - 2
            D4_ksi(i,j,1) = -eps_e*delta_t(i,j)*(1.0d0/metric_jacobian(i,j))*( Q_dis(i+2,j,1) &
                            - 4.0d0*Q_dis(i+1,j,1) - 4.0d0*Q_dis(i-1,j,1) + Q_dis(i-2,j,1) )
            D4_ksi(i,j,2) = -eps_e*delta_t(i,j)*(1.0d0/metric_jacobian(i,j))*( Q_dis(i+2,j,2) &
                            - 4.0d0*Q_dis(i+1,j,2) - 4.0d0*Q_dis(i-1,j,2) + Q_dis(i-2,j,2) )
            D4_ksi(i,j,3) = -eps_e*delta_t(i,j)*(1.0d0/metric_jacobian(i,j))*( Q_dis(i+2,j,3) &
                            - 4.0d0*Q_dis(i+1,j,3) - 4.0d0*Q_dis(i-1,j,3) + Q_dis(i-2,j,3) )
            D4_ksi(i,j,4) = -eps_e*delta_t(i,j)*(1.0d0/metric_jacobian(i,j))*( Q_dis(i+2,j,4) &
                            - 4.0d0*Q_dis(i+1,j,4) - 4.0d0*Q_dis(i-1,j,4) + Q_dis(i-2,j,4) )
            !
            !
            !         
            D4_eta(i,j,1) = -eps_e*delta_t(i,j)*(1.0d0/metric_jacobian(i,j))*( Q_dis(i,j+2,1) &
                            - 4.0d0*Q_dis(i,j+1,1) - 4.0d0*Q_dis(i,j-1,1) + Q_dis(i,j-2,1) )
            D4_eta(i,j,2) = -eps_e*delta_t(i,j)*(1.0d0/metric_jacobian(i,j))*( Q_dis(i,j+2,2) &
                            - 4.0d0*Q_dis(i,j+1,2) - 4.0d0*Q_dis(i,j-1,2) + Q_dis(i,j-2,2) )
            D4_eta(i,j,3) = -eps_e*delta_t(i,j)*(1.0d0/metric_jacobian(i,j))*( Q_dis(i,j+2,3) &
                            - 4.0d0*Q_dis(i,j+1,3) - 4.0d0*Q_dis(i,j-1,3) + Q_dis(i,j-2,3) )
            D4_eta(i,j,4) = -eps_e*delta_t(i,j)*(1.0d0/metric_jacobian(i,j))*( Q_dis(i,j+2,4) & 
                            - 4.0d0*Q_dis(i,j+1,4) - 4.0d0*Q_dis(i,j-1,4) + Q_dis(i,j-2,4) )
        end do
end do
!
! residue - pick the maximum residue of them 
! which means the norma_infinity
!
do j = 2, jmax - 1
        do i = 2, imax - 1
residue1(i,j) = 0.50d0*delta_ksi*(E_barra(i+1,j,1) - E_barra(i-1,j,1)) + D4_ksi(i,j,1) & 
                + 0.50d0*delta_eta*(F_barra(i,j+1,1) - F_barra(i,j-1,1)) + D4_eta(i,j,1)
        !
        !
residue2(i,j) = 0.50d0*delta_ksi*(E_barra(i+1,j,2) - E_barra(i-1,j,2)) + D4_ksi(i,j,2) &
                + 0.50d0*delta_eta*(F_barra(i,j+1,2) - F_barra(i,j-1,2)) + D4_eta(i,j,2)
        !
        !
residue3(i,j) = 0.50d0*delta_ksi*(E_barra(i+1,j,3) - E_barra(i-1,j,3)) + D4_ksi(i,j,3) & 
                + 0.50d0*delta_eta*(F_barra(i,j+1,3) - F_barra(i,j-1,3)) + D4_eta(i,j,3)
        !
        !
residue4(i,j) = 0.50d0*delta_ksi*(E_barra(i+1,j,4) - E_barra(i-1,j,4)) + D4_ksi(i,j,4) &
                + 0.50d0*delta_eta*(F_barra(i,j+1,4) - F_barra(i,j-1,4)) + D4_eta(i,j,4)
        !
        !
        max_residue = max(residue1(i,j), residue2(i,j), residue3(i,j), residue4(i,j) )
        end do
end do
!
! explicit time marching 
!
do j = 2, jmax - 1
        do i = 2, imax - 1
            Q_barra(i,j,1) = Q_barra(i,j,1) - delta_t(i,j)*residue1(i,j)
            Q_barra(i,j,2) = Q_barra(i,j,2) - delta_t(i,j)*residue2(i,j)
            Q_barra(i,j,3) = Q_barra(i,j,3) - delta_t(i,j)*residue3(i,j) 
            Q_barra(i,j,4) = Q_barra(i,j,4) - delta_t(i,j)*residue4(i,j) 
        end do
end do
!
!
!
deallocate(Q_dis, D4_ksi, D4_eta)
!
!
!
end subroutine euler_explicit
