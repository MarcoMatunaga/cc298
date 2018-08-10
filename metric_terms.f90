!
! metric terms
!
subroutine metric_terms
use vars
implicit none
!****************************************
! delta_eta e delta_ksi igual a 1: podemos colocar
! como delta_eta e delta_ksi como inputs
!****************************************
do j = 1, jmax
	do i = 1, imax 
		x_eta(i,j) = ( meshx(i,j+1) - meshx(i,j-1) )/( 2.0d0*delta_eta )
		y_eta(i,j) = ( meshy(i,j+1) - meshy(i,j-1) )/( 2.0d0*delta_eta )
		x_ksi(i,j) = ( meshx(i+1,j) - meshx(i-1,j) )/( 2.0d0*delta_ksi )
		y_ksi(i,j) = ( meshy(i+1,j) - meshy(i-1,j) )/( 2.0d0*delta_ksi )
		metric_jacobian = 1.0d0/(x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j))
		eta_x(i,j) = metric_jacobian*y_eta(i,j)
		eta_y(i,j) = -1.0d0*metric_jacobian*x_eta(i,j)
		ksi_x(i,j) = metric_jacobian*y_eta(i,j)
		ksi_y(i,j) = -1.0d0*metric_jacobian*y_ksi(i,j)
	end do
end do
!
!
!
end subroutine metric_terms
