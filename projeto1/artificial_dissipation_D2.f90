subroutine artificial_dissipation_D2
    use vars
    implicit none
!
!
eps_e = 5.00d0
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
! calculate in everybody 
!
do j = 2, jmax - 1
    do i = 2, imax - 1
            D4_ksi(i,j,1) = -eps_e*( Q_dis(i+1,j,1) - 2.0d0*Q_dis(i,j,1) &
                            + Q_dis(i-1,j,1) )

            D4_ksi(i,j,2) = -eps_e*( Q_dis(i+1,j,2) - 2.0d0*Q_dis(i,j,2) &
                            + Q_dis(i-1,j,2) )

            D4_ksi(i,j,3) = -eps_e*( Q_dis(i+1,j,3) - 2.0d0*Q_dis(i,j,3) &
                            + Q_dis(i-1,j,3) )

            D4_ksi(i,j,4) = -eps_e*( Q_dis(i+1,j,4) - 2.0d0*Q_dis(i,j,4) &
                            + Q_dis(i-1,j,4) )
            !
            !
            !         
            D4_eta(i,j,1) = -eps_e*( Q_dis(i,j+1,1) - 2.0d0*Q_dis(i,j,1) &
                            + Q_dis(i,j-1,1) )

            D4_eta(i,j,2) = -eps_e*( Q_dis(i,j+1,2) - 2.0d0*Q_dis(i,j,2) &
                            + Q_dis(i,j-1,2) )

            D4_eta(i,j,3) = -eps_e*( Q_dis(i,j+1,3) - 2.0d0*Q_dis(i,j,3) &
                            + Q_dis(i,j-1,3) )

            D4_eta(i,j,4) = -eps_e*( Q_dis(i,j+1,4) - 2.0d0*Q_dis(i,j,4) &
                            + Q_dis(i,j-1,4) )
    end do
end do
!
!
!
end subroutine artificial_dissipation_D2