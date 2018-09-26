module vars
implicit none
! flow, equations, variables vectors/matrixes
real(8),dimension(:,:),allocatable           :: T, p, u, v, q_vel, a
real(8),dimension(:,:),allocatable           :: U_contravariant, V_contravariant
real(8),dimension(:,:),allocatable           :: residue1, residue2, residue3, residue4 
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
integer(4)                                   :: iter, max_iter
! numerical variables, time marching
real(8),dimension(:,:),allocatable           :: delta_t
real(8)                                      :: delta_y, delta_x
real(8)                                      :: max_residue
real(8)                                      :: CFL
! mesh, metric terms, jacobians
real(8),dimension(:,:),allocatable           :: meshx, meshy
real(8),dimension(:,:),allocatable           :: eta_x, eta_y, ksi_x, ksi_y
real(8),dimension(:,:),allocatable           :: x_eta, y_eta, x_ksi, y_ksi
real(8),dimension(:,:),allocatable           :: metric_jacobian
real(8)                                      :: delta_eta, delta_ksi
!
integer(4)                                   :: nsave, total_sol
! artificial dissipation parameters
integer(4)                                   :: which_diss
real(8)                                      :: eps_dis_e, eps_dis_i
real(8),dimension(:,:,:), allocatable        :: Q_dis, D4_ksi, D4_eta

contains
!
!
!
    subroutine allocate_vars
        implicit none
!
!
!
        allocate(meshx(imax,jmax), meshy(imax,jmax))
        allocate(delta_t(imax,jmax))
        allocate(residue1(imax,jmax), residue2(imax,jmax))
        allocate(residue3(imax,jmax), residue4(imax,jmax))
        allocate(metric_jacobian(imax,jmax))
        allocate( U_contravariant(imax,jmax), V_contravariant(imax,jmax))
        allocate(p(imax,jmax), T(imax,jmax))
        allocate(u(imax,jmax), v(imax,jmax), a(imax,jmax), q_vel(imax,jmax))
        allocate(y_ksi(imax,jmax), y_eta(imax,jmax))
        allocate(x_ksi(imax,jmax), x_eta(imax,jmax))
        allocate(ksi_x(imax,jmax), ksi_y(imax,jmax))
        allocate(eta_x(imax,jmax), eta_y(imax,jmax))
        allocate(Q_barra(imax,jmax,dim), E_barra(imax,jmax,dim), F_barra(imax,jmax,dim))
        allocate(Q_dis(imax,jmax,dim), D4_ksi(imax,jmax,dim), D4_eta(imax,jmax,dim) )
!
!
!
    end subroutine allocate_vars
!
!
!
    subroutine deallocate_vars
        implicit none
        deallocate(meshx, meshy)
        deallocate(delta_t)
        deallocate(residue1, residue2)
        deallocate(residue3, residue4)
        deallocate(U_contravariant, V_contravariant)
        deallocate(p, T)
        deallocate(u, v, a, q_vel)
        deallocate(y_ksi, y_eta)
        deallocate(x_ksi, x_eta)
        deallocate(ksi_x, ksi_y)
        deallocate(eta_x, eta_y)
        deallocate(metric_jacobian)
        deallocate(Q_barra, E_barra, F_barra)
        deallocate(Q_dis, D4_ksi, D4_eta)
!
!
!
    end subroutine deallocate_vars
end module vars
