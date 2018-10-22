module diagonalization
    use vars
!
!
contains
!
!
    function diag_ksi(U_c,a_temp,ksi_x_temp,ksi_y_temp,size)
        implicit none
        integer(4),intent(in)                           :: size
        real(8),intent(in)                              :: U_c,a_temp,ksi_x_temp,ksi_y_temp
        real(8),dimension(size,size)                    :: diag_ksi
        diag_ksi = 0.0d0
        if (size == 4) then
            diag_ksi(1,1) = U_c
            diag_ksi(2,2) = U_c
            diag_ksi(3,3) = U_c + a_temp*sqrt(ksi_x_temp**2.0d0+ksi_y_temp**2.0d0) 
            diag_ksi(4,4) = U_c - a_temp*sqrt(ksi_x_temp**2.0d0+ksi_y_temp**2.0d0)
        endif
    end function diag_ksi
!
!
    function diag_eta(V_c,a_temp,eta_x_temp,eta_y_temp,size)
        implicit none
        integer(4),intent(in)                           :: size
        real(8),intent(in)                              :: V_c,a_temp,eta_x_temp,eta_y_temp
        real(8),dimension(size,size)                    :: diag_eta
        diag_eta = 0.0d0
        if (size == 4) then
            diag_eta(1,1) = V_c
            diag_eta(2,2) = V_c
            diag_eta(3,3) = V_c + a_temp*sqrt(eta_x_temp**2.0d0+eta_y_temp**2.0d0) 
            diag_eta(4,4) = V_c - a_temp*sqrt(eta_x_temp**2.0d0+eta_y_temp**2.0d0)
        endif
    end function diag_eta
!
!
    function T_ksi(u_temp,v_temp,rho_temp,a_temp,ksi_x_temp,ksi_y_temp,size)
        implicit none
        integer(4),intent(in)                            :: size
        real(8),intent(in)                               :: u_temp,v_temp,rho_temp,a_temp
        real(8),intent(in)                               :: ksi_x_temp, ksi_y_temp
        real(8),dimension(size,size)                     :: T_ksi              
        real(8)                                          :: alfa_t, phi_temp, theta_til
        real(8)                                          :: ksi_x_temp_til
        real(8)                                          :: ksi_y_temp_til  
        !
        !
        alfa_t         = rho_temp/(sqrt(2.0d0)*a_temp)
        ksi_x_temp_til = ksi_x_temp/(sqrt(ksi_x_temp**2.0d0+ksi_y_temp**2.0d0))
        ksi_y_temp_til = ksi_y_temp/(sqrt(ksi_x_temp**2.0d0+ksi_y_temp**2.0d0))
        phi_temp       = 0.50d0*(gama - 1.0d0)*(u_temp**2.0d0 + v_temp**2.0d0)
        theta_til      = ksi_x_temp_til*u_temp + ksi_y_temp_til*v_temp
        !
        !
        T_ksi(1,1) = 1.0d0
        T_ksi(1,2) = 0.0d0
        T_ksi(1,3) = alfa_t
        T_ksi(1,4) = alfa_t
        !
        !
        T_ksi(2,1) = u_temp
        T_ksi(2,2) = ksi_y_temp_til*rho_temp
        T_ksi(2,3) = alfa_t*(u_temp + ksi_x_temp_til*a_temp)
        T_ksi(2,4) = alfa_t*(u_temp - ksi_x_temp_til*a_temp)
        !
        !
        T_ksi(3,1) = v_temp
        T_ksi(3,2) = -ksi_x_temp_til*rho_temp
        T_ksi(3,3) = alfa_t*(v_temp + ksi_y_temp_til*a_temp)
        T_ksi(3,4) = alfa_t*(v_temp - ksi_y_temp_til*a_temp)
        !
        !
        T_ksi(4,1) = phi_temp**2.0d0/(gama - 1.0d0)
        T_ksi(4,2) = rho_temp*(ksi_y_temp_til*u_temp - ksi_x_temp_til*v_temp)
        T_ksi(4,3) = alfa_t*( (phi_temp**2.0d0 + a_temp**2.0d0)/(gama-1.0d0) + a_temp*theta_til )
        T_ksi(4,4) = alfa_t*( (phi_temp**2.0d0 + a_temp**2.0d0)/(gama-1.0d0) - a_temp*theta_til )
        !
        !
    end function T_ksi
!
!
    function inv_T_ksi(u_temp,v_temp,rho_temp,a_temp,ksi_x_temp,ksi_y_temp,size)
        implicit none
        integer(4),intent(in)                            :: size
        real(8),intent(in)                               :: u_temp,v_temp,rho_temp,a_temp
        real(8),intent(in)                               :: ksi_x_temp, ksi_y_temp
        real(8),dimension(size,size)                     :: inv_T_ksi              
        real(8)                                          :: beta_t, phi_temp, theta_til
        real(8)                                          :: ksi_x_temp_til
        real(8)                                          :: ksi_y_temp_til  
        !
        !
        beta_t         = 1.0d0/( sqrt(2.0d0)*a_temp*rho_temp )
        ksi_x_temp_til = ksi_x_temp/(sqrt(ksi_x_temp**2.0d0+ksi_y_temp**2.0d0))
        ksi_y_temp_til = ksi_y_temp/(sqrt(ksi_x_temp**2.0d0+ksi_y_temp**2.0d0))
        phi_temp       = 0.50d0*(gama - 1.0d0)*(u_temp**2.0d0 + v_temp**2.0d0)
        theta_til      = ksi_x_temp_til*u_temp + ksi_y_temp_til*v_temp
        !
        !
        inv_T_ksi(1,1) = ( 1.0d0 - phi_temp**2.0d0/a_temp**2.0d0)
        inv_T_ksi(1,2) = (gama - 1.0d0)*u_temp/a_temp**2.0d0
        inv_T_ksi(1,3) = (gama - 1.0d0)*v_temp/a_temp**2.0d0  
        inv_T_ksi(1,4) = -(gama - 1.0d0)/a_temp**2.0d0
        !
        !
        inv_T_ksi(2,1) = -(ksi_y_temp_til*u_temp - ksi_x_temp_til*v_temp)/rho_temp
        inv_T_ksi(2,2) = ksi_y_temp_til/rho_temp
        inv_T_ksi(2,3) = -ksi_x_temp_til/rho_temp 
        inv_T_ksi(2,4) = 0.0d0
        !
        !
        inv_T_ksi(3,1) = beta_t*(phi_temp**2.0d0 - a_temp*theta_til)
        inv_T_ksi(3,2) = beta_t*(ksi_x_temp_til*a_temp - (gama - 1.0d0)*u_temp )
        inv_T_ksi(3,3) = beta_t*(ksi_y_temp_til*a_temp - (gama - 1.0d0)*v_temp )
        inv_T_ksi(3,4) = beta_t*(gama - 1.0d0)
        !
        !
        inv_T_ksi(4,1) = beta_t*(phi_temp**2.0d0 + a_temp*theta_til)
        inv_T_ksi(4,2) = -beta_t*(ksi_x_temp_til*a_temp + gama - 1.0d0*u_temp )
        inv_T_ksi(4,3) = -beta_t*(ksi_y_temp_til*a_temp + gama - 1.0d0*v_temp )
        inv_T_ksi(4,4) = beta_t*(gama - 1.0d0)
        !
        !
    end function inv_T_ksi
!
!
    function T_eta(u_temp,v_temp,rho_temp,a_temp,eta_x_temp,eta_y_temp,size)
        implicit none
        integer(4),intent(in)                            :: size
        real(8),intent(in)                               :: u_temp,v_temp,rho_temp,a_temp
        real(8),intent(in)                               :: eta_x_temp, eta_y_temp
        real(8),dimension(size,size)                     :: T_eta
        real(8)                                          :: alfa_t, phi_temp, theta_til
        real(8)                                          :: eta_x_temp_til
        real(8)                                          :: eta_y_temp_til  
        !
        !
        alfa_t         = rho_temp/(sqrt(2.0d0)*a_temp)
        eta_x_temp_til = eta_x_temp/(sqrt(eta_x_temp**2.0d0+eta_y_temp**2.0d0))
        eta_y_temp_til = eta_y_temp/(sqrt(eta_x_temp**2.0d0+eta_y_temp**2.0d0))
        phi_temp       = 0.50d0*(gama - 1.0d0)*(u_temp**2.0d0 + v_temp**2.0d0)
        theta_til      = eta_x_temp_til*u_temp + eta_y_temp_til*v_temp
        !
        !
        T_eta(1,1) = 1.0d0
        T_eta(1,2) = 0.0d0
        T_eta(1,3) = alfa_t
        T_eta(1,4) = alfa_t
        !
        !
        T_eta(2,1) = u_temp
        T_eta(2,2) = eta_y_temp_til*rho_temp
        T_eta(2,3) = alfa_t*(u_temp + eta_x_temp_til*a_temp)
        T_eta(2,4) = alfa_t*(u_temp - eta_x_temp_til*a_temp)
        !
        !
        T_eta(3,1) = v_temp
        T_eta(3,2) = -eta_x_temp_til*rho_temp
        T_eta(3,3) = alfa_t*(v_temp + eta_y_temp_til*a_temp)
        T_eta(3,4) = alfa_t*(v_temp - eta_y_temp_til*a_temp)
        !
        !
        T_eta(4,1) = phi_temp**2.0d0/(gama - 1.0d0)
        T_eta(4,2) = rho_temp*(eta_y_temp_til*u_temp - eta_x_temp_til*v_temp)
        T_eta(4,3) = alfa_t*( (phi_temp**2.0d0 + a_temp**2.0d0)/(gama-1.0d0) + a_temp*theta_til )
        T_eta(4,4) = alfa_t*( (phi_temp**2.0d0 + a_temp**2.0d0)/(gama-1.0d0) - a_temp*theta_til )
        !
        !
    end function T_eta
!
!
    function inv_T_eta(u_temp,v_temp,rho_temp,a_temp,eta_x_temp,eta_y_temp,size)
        implicit none
        integer(4),intent(in)                            :: size
        real(8),intent(in)                               :: u_temp,v_temp,rho_temp,a_temp
        real(8),intent(in)                               :: eta_x_temp, eta_y_temp
        real(8),dimension(size,size)                     :: inv_T_eta
        real(8)                                          :: beta_t, phi_temp, theta_til
        real(8)                                          :: eta_x_temp_til
        real(8)                                          :: eta_y_temp_til  
        !
        !
        beta_t         = 1.0d0/( sqrt(2.0d0)*a_temp*rho_temp )
        eta_x_temp_til = eta_x_temp/(sqrt(eta_x_temp**2.0d0+eta_y_temp**2.0d0))
        eta_y_temp_til = eta_y_temp/(sqrt(eta_x_temp**2.0d0+eta_y_temp**2.0d0))
        phi_temp       = 0.50d0*(gama - 1.0d0)*(u_temp**2.0d0 + v_temp**2.0d0)
        theta_til      = eta_x_temp_til*u_temp + eta_y_temp_til*v_temp
        !
        !
        inv_T_eta(1,1) = ( 1.0d0 - phi_temp**2.0d0/a_temp**2.0d0)
        inv_T_eta(1,2) = (gama - 1.0d0)*u_temp/a_temp**2.0d0 
        inv_T_eta(1,3) = (gama - 1.0d0)*v_temp/a_temp**2.0d0  
        inv_T_eta(1,4) = -(gama - 1.0d0)/a_temp**2.0d0
        !
        !
        inv_T_eta(2,1) = -(eta_y_temp_til*u_temp - eta_x_temp_til*v_temp)/rho_temp 
        inv_T_eta(2,2) = eta_y_temp_til/rho_temp
        inv_T_eta(2,3) = -eta_x_temp_til/rho_temp  
        inv_T_eta(2,4) = 0.0d0
        !
        !
        inv_T_eta(3,1) = beta_t*(phi_temp**2.0d0 - a_temp*theta_til)
        inv_T_eta(3,2) = beta_t*(eta_x_temp_til*a_temp - (gama - 1.0d0)*u_temp )
        inv_T_eta(3,3) = beta_t*(eta_y_temp_til*a_temp - (gama - 1.0d0)*v_temp )
        inv_T_eta(3,4) = beta_t*(gama - 1.0d0)
        !
        !
        inv_T_eta(4,1) = beta_t*(phi_temp**2.0d0 + a_temp*theta_til)
        inv_T_eta(4,2) = -beta_t*(eta_x_temp_til*a_temp + gama - 1.0d0*u_temp )
        inv_T_eta(4,3) = -beta_t*(eta_y_temp_til*a_temp + gama - 1.0d0*v_temp )
        inv_T_eta(4,4) = beta_t*(gama - 1.0d0)
        !
        !
    end function inv_T_eta
!
!
    function N_matrix(ksi_x_temp,ksi_y_temp,eta_x_temp,eta_y_temp,size)
        implicit none
        integer(4),intent(in)                            :: size
        real(8),intent(in)                               :: eta_x_temp,eta_y_temp
        real(8),intent(in)                               :: ksi_x_temp,ksi_y_temp
        real(8)                                          :: ksi_x_temp_til
        real(8)                                          :: ksi_y_temp_til
        real(8)                                          :: eta_x_temp_til
        real(8)                                          :: eta_y_temp_til
        real(8)                                          :: m1,m2
        real(8),dimension(size,size)                     :: N_matrix
        !
        !
        ksi_x_temp_til = ksi_x_temp/(sqrt(ksi_x_temp**2.0d0+ksi_y_temp**2.0d0))
        ksi_y_temp_til = ksi_y_temp/(sqrt(ksi_x_temp**2.0d0+ksi_y_temp**2.0d0))
        eta_x_temp_til = eta_x_temp/(sqrt(eta_x_temp**2.0d0+eta_y_temp**2.0d0))
        eta_y_temp_til = eta_y_temp/(sqrt(eta_x_temp**2.0d0+eta_y_temp**2.0d0))
        m1             = (ksi_x_temp_til*eta_x_temp_til + ksi_y_temp_til*eta_y_temp_til)
        m2             = (ksi_x_temp_til*eta_x_temp_til - ksi_y_temp_til*eta_y_temp_til)
        !
        !
        N_matrix(1,1) = 1.0d0
        N_matrix(1,2) = 0.0d0
        N_matrix(1,3) = 0.0d0
        N_matrix(1,4) = 0.0d0
        !
        !
        N_matrix(2,1) = 0.0d0
        N_matrix(2,2) = m1
        N_matrix(2,3) = -(1.0d0/sqrt(2.0d0))*m2
        N_matrix(2,4) = (1.0d0/sqrt(2.0d0))*m2
        !
        !
        N_matrix(3,1) = 0.0d0
        N_matrix(3,2) = (1.0d0/sqrt(2.0d0))*m2
        N_matrix(3,3) = (1.0d0/sqrt(2.0d0))**2.0d0*(1.0d0 + m1)
        N_matrix(3,4) = (1.0d0/sqrt(2.0d0))**2.0d0*(1.0d0 - m1)
        !
        !
        N_matrix(4,1) = 0.0d0
        N_matrix(4,2) = -(1.0d0/sqrt(2.0d0))*m2
        N_matrix(4,3) = (1.0d0/sqrt(2.0d0))**2.0d0*(1.0d0 - m1)
        N_matrix(4,4) = (1.0d0/sqrt(2.0d0))**2.0d0*(1.0d0 + m1)
        !
        !
    end function N_matrix
!
!
    function inv_N_matrix(ksi_x_temp,ksi_y_temp,eta_x_temp,eta_y_temp,size)
        implicit none
        integer(4),intent(in)                            :: size
        real(8),intent(in)                               :: eta_x_temp,eta_y_temp
        real(8),intent(in)                               :: ksi_x_temp,ksi_y_temp
        real(8)                                          :: ksi_x_temp_til
        real(8)                                          :: ksi_y_temp_til
        real(8)                                          :: eta_x_temp_til
        real(8)                                          :: eta_y_temp_til
        real(8)                                          :: m1,m2
        real(8),dimension(size,size)                     :: inv_N_matrix
        !
        !
        ksi_x_temp_til = ksi_x_temp/(sqrt(ksi_x_temp**2.0d0+ksi_y_temp**2.0d0))
        ksi_y_temp_til = ksi_y_temp/(sqrt(ksi_x_temp**2.0d0+ksi_y_temp**2.0d0))
        eta_x_temp_til = eta_x_temp/(sqrt(eta_x_temp**2.0d0+eta_y_temp**2.0d0))
        eta_y_temp_til = eta_y_temp/(sqrt(eta_x_temp**2.0d0+eta_y_temp**2.0d0))
        m1             = (ksi_x_temp_til*eta_x_temp_til + ksi_y_temp_til*eta_y_temp_til)
        m2             = (ksi_x_temp_til*eta_y_temp_til - ksi_y_temp_til*eta_x_temp_til)
        !
        !
        inv_N_matrix(1,1) = 1.0d0
        inv_N_matrix(1,2) = 0.0d0
        inv_N_matrix(1,3) = 0.0d0
        inv_N_matrix(1,4) = 0.0d0
        !
        !
        inv_N_matrix(2,1) = 0.0d0
        inv_N_matrix(2,2) = m1
        inv_N_matrix(2,3) = (1.0d0/sqrt(2.0d0))*m2
        inv_N_matrix(2,4) = -(1.0d0/sqrt(2.0d0))*m2
        !
        !
        inv_N_matrix(3,1) = 0.0d0
        inv_N_matrix(3,2) = -(1.0d0/sqrt(2.0d0))*m2
        inv_N_matrix(3,3) = (1.0d0/sqrt(2.0d0))**2.0d0*(1.0d0 + m1)
        inv_N_matrix(3,4) = (1.0d0/sqrt(2.0d0))**2.0d0*(1.0d0 - m1)
        !
        !
        inv_N_matrix(4,1) = 0.0d0
        inv_N_matrix(4,2) = (1.0d0/sqrt(2.0d0))*m2
        inv_N_matrix(4,3) = (1.0d0/sqrt(2.0d0))**2.0d0*(1.0d0 - m1)
        inv_N_matrix(4,4) = (1.0d0/sqrt(2.0d0))**2.0d0*(1.0d0 + m1)
        !
        !
    end function inv_N_matrix
!
!
end module diagonalization
