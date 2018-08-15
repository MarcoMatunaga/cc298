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
do j = 2, jmax - 1
	do i = 2, imax - 1 
		x_eta(i,j) = ( meshx(i,j+1) - meshx(i,j-1) )/( 2.0d0*delta_eta )
		y_eta(i,j) = ( meshy(i,j+1) - meshy(i,j-1) )/( 2.0d0*delta_eta )
		x_ksi(i,j) = ( meshx(i+1,j) - meshx(i-1,j) )/( 2.0d0*delta_ksi )
		y_ksi(i,j) = ( meshy(i+1,j) - meshy(i-1,j) )/( 2.0d0*delta_ksi )
	end do
end do
!
!
i = 1
do j = 1, jmax
	x_ksi(i,j) = ( meshx(i+1,j) - meshx(i,j) )/( delta_ksi )
	y_ksi(i,j) = ( meshy(i+1,j) - meshy(i,j) )/( delta_ksi )
end do
!
!
i = imax
do j = 1, jmax
	x_ksi(i,j) = ( meshx(i,j) - meshx(i-1,j) )/( delta_ksi )
	y_ksi(i,j) = ( meshy(i,j) - meshy(i-1,j) )/( delta_ksi )
end do
!
!
j = 1
do i = 1, imax
	x_eta(i,j) = ( meshx(i,j+1) - meshx(i,j) )/( delta_eta )
	y_eta(i,j) = ( meshy(i,j+1) - meshy(i,j) )/( delta_eta )
end do
!
!
j = jmax
do i = 1, imax
	x_eta(i,j) = ( meshx(i,j) - meshx(i,j-1) )/( delta_eta )
	y_eta(i,j) = ( meshy(i,j) - meshy(i,j-1) )/( delta_eta )		
end do
!
!
!
do j = 1, jmax
	do i = 1, imax
		metric_jacobian(i,j) = 1.0d0/(x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j))
		eta_x(i,j) = -1.0d0*metric_jacobian(i,j)*y_ksi(i,j)
		eta_y(i,j) = metric_jacobian(i,j)*x_ksi(i,j)
		ksi_x(i,j) = metric_jacobian(i,j)*y_eta(i,j)
		ksi_y(i,j) = -1.0d0*metric_jacobian(i,j)*x_eta(i,j)
	end do
end do
!
!
!
end subroutine metric_terms
