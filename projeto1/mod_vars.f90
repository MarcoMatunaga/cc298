module vars
implicit none

! flow, equations, variables vectors/matrixes
real(8),dimension(:,:),allocatable           :: U_contravariant, V_contravariant
real(8),dimension(:,:,:),allocatable         :: residue
real(8),dimension(:,:,:),allocatable         :: Q_barra, E_barra, F_barra 
real(8),dimension(:,:,:),allocatable         :: E_pos, E_neg, F_pos, F_neg
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
real(8)                                      :: max_residue, res_conv
real(8)                                      :: CFL
integer(4)                                   :: ramp, CFL_ramp_size


! mesh, metric terms, jacobians
real(8),dimension(:,:),allocatable           :: meshx, meshy
real(8),dimension(:,:),allocatable           :: eta_x, eta_y, ksi_x, ksi_y
real(8),dimension(:,:),allocatable           :: x_eta, y_eta, x_ksi, y_ksi
real(8),dimension(:,:),allocatable           :: metric_jacobian
real(8)                                      :: delta_eta, delta_ksi

integer(4)                                   :: nsave, total_sol

! time marching method
integer(4)                                   :: time_method
integer(4)                                   :: which_boundary
    
! artificial dissipation parameters
integer(4)                                   :: which_diss
real(8)                                      :: eps_dis_e, eps_dis_i, dis_factor, k2, k4
real(8)                                      :: Dnonlinear_plus, Dnonlinear_minus
real(8),dimension(:,:,:),allocatable         :: Q_dis, D4_ksi, D4_eta

! non linear dissipation parameters
real(8),dimension(:,:),allocatable           :: nu_dis_eta, nu_dis_ksi
real(8),dimension(:,:),allocatable           :: eps2_ksi, eps4_ksi
real(8),dimension(:,:),allocatable           :: eps2_eta, eps4_eta
real(8),dimension(:,:),allocatable           :: sigma_ksi, sigma_eta
    
    contains

    subroutine allocate_vars
        implicit none
        allocate(meshx(imax,jmax), meshy(imax,jmax))
        allocate(delta_t(imax,jmax))
        allocate(residue(imax,jmax,dim))
        allocate(metric_jacobian(imax,jmax))
        allocate(U_contravariant(imax,jmax), V_contravariant(imax,jmax))
        allocate(y_ksi(imax,jmax), y_eta(imax,jmax))
        allocate(x_ksi(imax,jmax), x_eta(imax,jmax))
        allocate(ksi_x(imax,jmax), ksi_y(imax,jmax))
        allocate(eta_x(imax,jmax), eta_y(imax,jmax))
        allocate(Q_barra(imax,jmax,dim), E_barra(imax,jmax,dim), F_barra(imax,jmax,dim))
        allocate(Q_dis(imax,jmax,dim), D4_ksi(imax,jmax,dim), D4_eta(imax,jmax,dim) )
        allocate(E_pos(imax,jmax,dim),E_neg(imax,jmax,dim))
        allocate(F_pos(imax,jmax,dim),F_neg(imax,jmax,dim))
    end subroutine allocate_vars

    subroutine deallocate_vars
        implicit none
        deallocate(meshx, meshy)
        deallocate(delta_t)
        deallocate(residue)
        deallocate(U_contravariant, V_contravariant)
        deallocate(y_ksi, y_eta)
        deallocate(x_ksi, x_eta)
        deallocate(ksi_x, ksi_y)
        deallocate(eta_x, eta_y)
        deallocate(metric_jacobian)
        deallocate(Q_barra, E_barra, F_barra)
        deallocate(Q_dis, D4_ksi, D4_eta)
        deallocate(E_pos,E_neg)
        deallocate(F_pos,F_neg)
    end subroutine deallocate_vars

    subroutine allocate_vars_non_linear
        implicit none
        allocate(nu_dis_ksi(imax,jmax),nu_dis_eta(imax,jmax))
        allocate(eps2_ksi(imax,jmax),eps4_ksi(imax,jmax))
        allocate(eps2_eta(imax,jmax),eps4_eta(imax,jmax))
        allocate(sigma_ksi(imax,jmax), sigma_eta(imax,jmax))
    end subroutine allocate_vars_non_linear

    subroutine deallocate_vars_non_linear
        implicit none
        deallocate(nu_dis_ksi,nu_dis_eta)
        deallocate(eps2_ksi,eps4_ksi)
        deallocate(eps2_eta,eps4_eta)
        deallocate(sigma_ksi,sigma_eta)
    end subroutine deallocate_vars_non_linear

end module vars
