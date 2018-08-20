!
! euler explicit
!
subroutine euler_explicit
use vars
implicit none
!
!
!
do j = 2, jmax - 1
        do i = 2, imax - 1
residue1 = 0.50d0*delta_ksi*(E_barra(i+1,j,1) - E_barra(i-1,j,1)) + 0.50d0*delta_eta*(F_barra(i,j+1,1) - F_barra(i,j-1,1))
residue2 = 0.50d0*delta_ksi*(E_barra(i+1,j,2) - E_barra(i-1,j,2)) + 0.50d0*delta_eta*(F_barra(i,j+1,2) - F_barra(i,j-1,2))
residue3 = 0.50d0*delta_ksi*(E_barra(i+1,j,3) - E_barra(i-1,j,3)) + 0.50d0*delta_eta*(F_barra(i,j+1,3) - F_barra(i,j-1,3)) 
residue4 = 0.50d0*delta_ksi*(E_barra(i+1,j,4) - E_barra(i-1,j,4)) + 0.50d0*delta_eta*(F_barra(i,j+1,4) - F_barra(i,j-1,4))
        end do
end do
!
!
end subroutine euler_explicit