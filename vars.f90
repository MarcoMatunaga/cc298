module vars
implicit none
! flow, equations, variables vectors/matrixes
real(8),dimension(:,:),allocatable           :: T, p, u, v, q_vel, a
real(8),dimension(:,:,:),allocatable         :: Q, E, F
real(8),dimension(:,:,:),allocatable         :: Q_barra, E_barra, F_barra 
real(8)                                      :: T_total, p_total
! vector dimension 
integer(4)                                   :: dim 
! flow properties
real(8)                                      :: gama, c_v, R, a_cr
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

contains
    subroutine allocate_vars
        implicit none
!
!
!
        allocate(meshx(imax,jmax), meshy(imax,jmax))
        allocate(Q(imax,jmax,dim))
        allocate(phi_jacobian(imax,jmax), theta_jacobian(imax,jmax))
        allocate(a1_jacobian(imax,jmax), metric_jacobian(imax,jmax))
        allocate(p(imax,jmax), T(imax,jmax))
        allocate(u(imax,jmax), v(imax,jmax), a(imax,jmax), q_vel(imax,jmax))
        allocate(y_ksi(imax,jmax), y_eta(imax,jmax))
        allocate(x_ksi(imax,jmax), x_eta(imax,jmax))
        allocate(ksi_x(imax,jmax), ksi_y(imax,jmax))
        allocate(eta_x(imax,jmax), eta_y(imax,jmax))
        allocate(A_barra(imax,jmax), B_barra(imax,jmax), M_barra(imax,jmax))
        allocate(E_barra(imax,jmax,dim), F_barra(imax,jmax,dim))
!
!
!
    end subroutine allocate_vars
!
!
!
    subroutine deallocate_vars
        implicit none
        deallocate(meshx, meshy, Q)
        deallocate(p, T)
        deallocate(u, v, a, q_vel)
        deallocate(y_ksi, y_eta)
        deallocate(x_ksi, x_eta)
        deallocate(ksi_x, ksi_y)
        deallocate(eta_x, eta_y)
        deallocate(phi_jacobian, theta_jacobian, a1_jacobian)
        deallocate(metric_jacobian)
        deallocate(A_barra, B_barra, M_barra)
        deallocate(E_barra, F_barra)
!
!
!
end subroutine deallocate_vars
end module vars