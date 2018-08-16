module vars
implicit none
! flow, equations, variables vectors/matrixes
real(8),dimension(:,:),allocatable           :: T, p, u, v, q_vel, a, a_cr
real(8),dimension(:,:,:),allocatable         :: Q, E, F
real(8),dimension(:,:,:),allocatable         :: Q_barra, E_barra, F_barra 
real(8)                                      :: T_total, p_total
! vector dimension 
integer(4)                                   :: dim 
! flow properties
real(8)                                      :: gama, c_v, R
! flow constants
real(8)                                      :: theta 
! math constants
real(8)                                      :: pi, dummy
! indices dos vetores --> geometria, malha
integer(4)                                   :: i,j,k
integer(4)                                   :: imax,jmax,kmax
! numerical variables, time marching
real(8)                                      :: delta_t, delta_y, delta_x
real(8)                                      :: residue1, residue2, residue3, residue4
integer(4)                                   :: CFL
! mesh, metric terms, jacobians
real(8),dimension(:,:),allocatable           :: meshx, meshy
real(8),dimension(:,:),allocatable           :: eta_x, eta_y, ksi_x, ksi_y
real(8),dimension(:,:),allocatable           :: x_eta, y_eta, x_ksi, y_ksi
real(8),dimension(:,:),allocatable           :: metric_jacobian
real(8),dimension(:,:),allocatable           :: phi_jacobian, theta_jacobian, a1_jacobian
real(8),dimension(:,:),allocatable           :: A_barra, B_barra, M_barra
real(8)                                      :: delta_eta, delta_ksi


    dim = 4
    theta = 0.0d0*(180.0d0/pi)
    gama = 1.4d0
    pi = 3.1415d0
    delta_eta = 1.0d0
    delta_ksi = 1.0d0
    !
    ! here we consider the value of cv as (5/2)*R
    ! R is the universal perfect gas constant
    !
    R = 287.053d0
    c_v = 2.5d0*R
    !T_total = 0.555556d0*531.2d0
    T_total = 294.8d0
    !p_total = 47.880258888889d0*2117.0d0
    p_total = 101360.0d0

end module vars
