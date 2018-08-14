module vars
implicit none
! flow, equations, variables vectors/matrixes
real(8),dimension(:,:),allocatable  		 :: T, p, u, v, rho, q_vel, a, a_cr
real(8),dimension(:,:,:),allocatable 		 :: Q 
real(8)						 :: T_total, p_total
! flow properties
real(8)						 :: gama, c_v, R
! flow constants
real(8)						 :: theta 
! math constants
real(8)						 :: pi, dummy
! indices dos vetores --> geometria, malha
real(8),dimension(:,:),allocatable               :: meshx, meshy
integer(4)				         :: i,j,k
integer(4)					 :: imax,jmax,kmax
! numerical variables
real(8)						 :: delta_t, delta_y, delta_x
real(8)						 :: residue1, residue2, residue3, residue4
integer(4)					 :: CFL
! metric terms
real(8),dimension(:,:),allocatable               :: eta_x, eta_y, ksi_x, ksi_y
real(8),dimension(:,:),allocatable               :: x_eta, y_eta, x_ksi, y_ksi
real(8)						 :: metric_jacobian
real(8)						 :: delta_eta, delta_ksi
end module vars
