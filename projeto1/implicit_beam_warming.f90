!
! implicit_beam_warming
!
subroutine implicit_beam_warming
    use vars
    implicit none
!
real(8),dimension(:,:),allocatable           :: A_barra, B_barra, M_barra
!
!
!
allocate(A_barra(dim,dim), B_barra(dim,dim), M_barra(dim,dim))
!
! residue - pick the maximum residue of them 
! which means the norma_infinity
!
max_residue = -1.0d0
do j = 2, jmax - 1
        do i = 2, imax - 1
!
! remember that the Q_barra is multiplied by metric_jacobian
! when i call the jacobian i divide energy by rho which 
! cancels the metric_jacobian
!
call jacobian(u(i,j),v(i,j),Q_barra(i,j,4),Q_barra(i,j,1),ksi_x(i,j),ksi_y(i,j),eta_x(i,j),eta_y(i,j),dim,A_barra,B_barra)
residue1(i,j) = 0.50d0*(E_barra(i+1,j,1) - E_barra(i-1,j,1))  & 
              + 0.50d0*(F_barra(i,j+1,1) - F_barra(i,j-1,1)) 
        !
        !
residue2(i,j) = 0.50d0*(E_barra(i+1,j,2) - E_barra(i-1,j,2))  &
              + 0.50d0*(F_barra(i,j+1,2) - F_barra(i,j-1,2)) 
        !
        !
residue3(i,j) = 0.50d0*(E_barra(i+1,j,3) - E_barra(i-1,j,3)) & 
              + 0.50d0*(F_barra(i,j+1,3) - F_barra(i,j-1,3)) 
        !
        !
residue4(i,j) = 0.50d0*(E_barra(i+1,j,4) - E_barra(i-1,j,4)) &
              + 0.50d0*(F_barra(i,j+1,4) - F_barra(i,j-1,4)) 
        !
        !
        residue1(i,j) = residue1(i,j) + metric_jacobian(i,j)*(D4_ksi(i,j,1) + D4_eta(i,j,1))
        residue2(i,j) = residue2(i,j) + metric_jacobian(i,j)*(D4_ksi(i,j,2) + D4_eta(i,j,2))
        residue3(i,j) = residue3(i,j) + metric_jacobian(i,j)*(D4_ksi(i,j,3) + D4_eta(i,j,3))
        residue4(i,j) = residue4(i,j) + metric_jacobian(i,j)*(D4_ksi(i,j,4) + D4_eta(i,j,4))
        !
        !
        residue1(i,j) = delta_t(i,j)*residue1(i,j)
        residue2(i,j) = delta_t(i,j)*residue2(i,j)
        residue3(i,j) = delta_t(i,j)*residue3(i,j)
        residue4(i,j) = delta_t(i,j)*residue4(i,j)
        !
        !
        !
        if ( abs(residue1(i,j)) > max_residue ) max_residue = log10(abs(residue1(i,j)))
        if ( abs(residue2(i,j)) > max_residue ) max_residue = log10(abs(residue2(i,j)))
        if ( abs(residue3(i,j)) > max_residue ) max_residue = log10(abs(residue3(i,j)))
        if ( abs(residue4(i,j)) > max_residue ) max_residue = log10(abs(residue4(i,j)))
        !
        !
        end do
end do
!
! explicit time marching 
!
do j = 2, jmax - 1
        do i = 2, imax - 1
            Q_barra(i,j,1) = Q_barra(i,j,1) - residue1(i,j)
            Q_barra(i,j,2) = Q_barra(i,j,2) - residue2(i,j)
            Q_barra(i,j,3) = Q_barra(i,j,3) - residue3(i,j) 
            Q_barra(i,j,4) = Q_barra(i,j,4) - residue4(i,j) 
        end do
end do
!
!
!
deallocate(A_barra, B_barra, M_barra)
!
!
!
end subroutine implicit_beam_warming
