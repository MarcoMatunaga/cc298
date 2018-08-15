module vars
implicit none
! flow, equations, variables vectors/matrixes
real(8),dimension(:,:),allocatable  		 :: T, p, u, v, rho, q_vel, a, a_cr
real(8),dimension(:,:,:),allocatable 		 :: Q, E, F
real(8),dimension(:,:,:),allocatable 		 :: Q_barra, E_barra, F_barra 
real(8)						 :: T_total, p_total
! flow properties
real(8)						 :: gama, c_v, R
! flow constants
real(8)						 :: theta 
! math constants
real(8)						 :: pi, dummy
! indices dos vetores --> geometria, malha
integer(4)				         :: i,j,k
integer(4)					 :: imax,jmax,kmax
! numerical variables, time marching
real(8)						 :: delta_t, delta_y, delta_x
real(8)						 :: residue1, residue2, residue3, residue4
integer(4)					 :: CFL
! mesh, metric terms, jacobians
real(8),dimension(:,:),allocatable               :: meshx, meshy
real(8),dimension(:,:),allocatable               :: eta_x, eta_y, ksi_x, ksi_y
real(8),dimension(:,:),allocatable               :: x_eta, y_eta, x_ksi, y_ksi
real(8),dimension(:,:),allocatable		 :: metric_jacobian
real(8),dimension(:,:),allocatable               :: phi_jacobian, theta_jacobian, a1_jacobian
real(8),dimension(:,:),allocatable		 :: A_barra, B_barra, M_barra
real(8)						 :: delta_eta, delta_ksi
end module vars
