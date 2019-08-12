module fvs_van_leer
    use vars
    implicit none
    real :: 
    
    contains
    
    ! subroutine eigenvalues_vl(args)
    !     implicit none
    !     real :: args
        
    ! end subroutine eigenvalues_vl

    subroutine calculate_fluxes_vl(flux_ksi_pos,flux_ksi_neg,flux_eta_pos,flux_eta_neg)
        implicit none
        real(8),dimension(imax,jmax,dim),intent(out)     :: flux_ksi_pos, flux_ksi_neg
        real(8),dimension(imax,jmax,dim),intent(out)     :: flux_eta_pos, flux_eta_neg  
        integer(4)                                       :: rho,u,v,p,a
        real(8)                                          :: fmass_pos, fmass_neg
        real(8)                                          :: const

        do aux_j = 1, jmax
            do aux_i = 1, imax
        ! the variables are initialized in the module vars_vl
                rho(aux_i,aux_j) = Q_barra(aux_i,aux_j,1)/metric_jacobian(aux_i,aux_j)
                u(aux_i,aux_j) = Q_barra(aux_i,aux_j,2)/Q_barra(aux_i,aux_j,1)
                v(aux_i,aux_j) = Q_barra(aux_i,aux_j,3)/Q_barra(aux_i,aux_j,1)
                p(aux_i,aux_j) = (gama-1.0d0) * (Q_barra(aux_i,aux_j,4)/metric_jacobian(aux_i,aux_j) & 
                                 - 0.5d0*rho(aux_i,aux_j)*( (u(aux_i,aux_j))**2.0d0 + (v(aux_i,aux_j))**2.0d0))
                a(aux_i,aux_j) = sqrt(gama*p(aux_i,aux_j)/rho(aux_i,aux_j))  
        ! compute Epos
                if ( u(aux_i,aux_j)/a(aux_i,aux_j) >= 1.0d0 ) then
                    call fluxes_curvilinear_vl(u(aux_i,aux_j),v(aux_i,aux_j),p(aux_i,aux_j),aux_i,aux_j)
                    flux_ksi_pos = E_barra
                    flux_ksi_neg = 0.0d0
                else if (u(aux_i,aux_j)/a(aux_i,aux_j) < 1.0d0) then 
                    call fluxes_curvilinear_vl(u(aux_i,aux_j),v(aux_i,aux_j),p(aux_i,aux_j),aux_i,aux_j)
                    flux_ksi_pos = 0.0d0
                    flux_ksi_neg = E_barra
                else if (-1.0d0 < u(aux_i,aux_j)/a(aux_i,aux_j) < 1.0e-5)
                    const = metric_jacobian(aux_i,aux_j)*rho(aux_i,aux_j)*a(aux_i,aux_j)/(2.0d0*gama)  
                    fb = (gama-1.0d0)*u+2.0d0*a(aux_i,aux_j)   
                    flux_ksi_pos(aux_i,aux_j,1) = fmass_pos
                    flux_ksi_pos(aux_i,aux_j,2) = fmass_pos*(1.0d0/gama)*((gama-1.0d0)*u+2.0d0*a(aux_i,aux_j))
                    flux_ksi_pos(aux_i,aux_j,3) = fmass_pos*v(aux_i,aux_j)
                    flux_ksi_pos(aux_i,aux_j,4) = fmass_pos*(fb**2.0d0/(2.0d0*(gama**2.0d0 - 1.0d0)) + 0.50d0*v(aux_i,aux_j)**2.0d0)
                    
                    flux_ksi_neg(aux_i,aux_j,1) = 
                    flux_ksi_neg(aux_i,aux_j,2) = 
                    flux_ksi_neg(aux_i,aux_j,3) = 
                    flux_ksi_neg(aux_i,aux_j,4) = 
                end if
            end do 
        end do

    end subroutine calculate_fluxes_vl

    subroutine fluxes_curvilinear_vl(u_flux,v_flux,p_flux,ii,jj)
        implicit none
        real(8),intent(in)                       :: u_flux, v_flux, p_flux
        integer(4),intent(in)                    :: ii, jj
        !********
        !modify: calculate only in the inner points 
        ! e pass the boundary points to the boundary subroutine
        !********
        !********
        !change the way you calculate
        !********
        E_barra = 0.0d0
        F_barra = 0.0d0
        
        if ( jj /= 1 .and. jj/= jmax .and. ii /= 1 .and. ii /= imax ) then
            U_contravariant(ii,jj) = u_flux*ksi_x(ii,jj) + v_flux*ksi_y(ii,jj)
            V_contravariant(ii,jj) = u_flux*eta_x(ii,jj) + v_flux*eta_y(ii,jj)
        endif

        !  fluxes in curvilinear coordinates - ksi direction
        E_barra(ii,jj,1) = Q_barra(ii,jj,1)*U_contravariant(ii,jj) 
        E_barra(ii,jj,2) = Q_barra(ii,jj,2)*U_contravariant(ii,jj) + metric_jacobian(ii,jj)*p_flux*ksi_x(ii,jj)
        E_barra(ii,jj,3) = Q_barra(ii,jj,3)*U_contravariant(ii,jj) + metric_jacobian(ii,jj)*p_flux*ksi_y(ii,jj)
        E_barra(ii,jj,4) = metric_jacobian(ii,jj)*( (Q_barra(ii,jj,4)/metric_jacobian(ii,jj)) + p_flux )*U_contravariant(ii,jj)
        
        !  fluxes in curvilinear coordinates - eta direction
        F_barra(ii,jj,1) = Q_barra(ii,jj,1)*V_contravariant(ii,jj) 
        F_barra(ii,jj,2) = Q_barra(ii,jj,2)*V_contravariant(ii,jj) + metric_jacobian(ii,jj)*p_flux*eta_x(ii,jj)
        F_barra(ii,jj,3) = Q_barra(ii,jj,3)*V_contravariant(ii,jj) + metric_jacobian(ii,jj)*p_flux*eta_y(ii,jj)
        F_barra(ii,jj,4) = metric_jacobian(ii,jj)*( (Q_barra(ii,jj,4)/metric_jacobian(ii,jj)) + p_flux )*V_contravariant(ii,jj)

    end subroutine fluxes_curvilinear_vl

end module fvs_van_leer