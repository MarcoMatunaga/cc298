!
! compute the jacobians
!
subroutine jacobian_ksi(u_jac,v_jac,e_jac,rho_jac,ksi_x_jac,ksi_y_jac,size,A_jac)
    use vars
    implicit none
    real(8), intent(in)                              :: u_jac,v_jac,e_jac,rho_jac,ksi_x_jac,ksi_y_jac    
    integer(4), intent(in)                           :: size
    real(8),dimension(size,size),intent(out)         :: A_jac
    real(8)                                          :: phi_jacobian, theta_jacobian_ksi, a1_jacobian
!
! calculate some terms present in the Euler Jacobians
!
phi_jacobian       = sqrt( 0.5d0*(gama - 1.0d0)*(u_jac**2.0d0 + v_jac**2.0d0) )
a1_jacobian        = gama*( e_jac/rho_jac ) - phi_jacobian**2.0d0
theta_jacobian_ksi = ksi_x_jac*u_jac + ksi_y_jac*v_jac
!
! Euler jacobian in the ksi-direction
!
A_jac(1,1) = 0.0d0 
A_jac(1,2) = ksi_x_jac 
A_jac(1,3) = ksi_y_jac
A_jac(1,4) = 0.0d0
!
A_jac(2,1) = -u_jac*theta_jacobian_ksi + ksi_x_jac*phi_jacobian**2.0d0  
A_jac(2,2) = theta_jacobian_ksi - (gama - 2.0d0)*ksi_x_jac*u_jac
A_jac(2,3) = ksi_y_jac*u_jac - (gama - 1.0d0)*ksi_x_jac*v_jac
A_jac(2,4) = (gama - 1.0d0)*ksi_x_jac
!
A_jac(3,1) = -v_jac*theta_jacobian_ksi + ksi_y_jac*phi_jacobian**2.0d0
A_jac(3,2) = ksi_x_jac*v_jac - (gama - 1.0d0)*ksi_y_jac*u_jac
A_jac(3,3) = theta_jacobian_ksi - (gama - 2.0d0)*ksi_y_jac*v_jac
A_jac(3,4) = (gama - 1.0d0)*ksi_y_jac
!
A_jac(4,1) = theta_jacobian_ksi*(phi_jacobian**2.0d0 - a1_jacobian)
A_jac(4,2) = ksi_x_jac*a1_jacobian - (gama - 1.0d0)*u_jac*theta_jacobian_ksi
A_jac(4,3) = ksi_y_jac*a1_jacobian - (gama - 1.0d0)*v_jac*theta_jacobian_ksi
A_jac(4,4) = gama*theta_jacobian_ksi
!
!
end subroutine jacobian_ksi
