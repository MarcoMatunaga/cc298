module fluxes_pos_neg
    use vars
    use vars_sw
    use diagonalization
    
    contains
    
    function fluxes_split(k1,k2,eigen1,eigen3,eigen4,ind_i,ind_j)
        implicit none
        real(8),intent(in)                  :: k1,k2
        real(8),intent(in)                  :: eigen1,eigen3,eigen4
        integer(4),intent(in)               :: ind_i,ind_j
        real(8),dimension(dim)              :: fluxes_split
        real(8)                             :: temp
        real(8)                             :: k1_til, k2_til, wII

        k1_til = k1/sqrt(k1**2.0d0 + k2**2.0d0)
        k2_til = k2/sqrt(k1**2.0d0 + k2**2.0d0)
        
        fluxes_split(1) = 2.0d0*(gama - 1.0d0)*eigen1 + eigen3 + eigen4

        fluxes_split(2) = 2.0d0*(gama - 1.0d0)*eigen1*u(ind_i,ind_j) + &
                                eigen3*(u(ind_i,ind_j) + a(ind_i,ind_j)*k1_til) + &
                                eigen4*(u(ind_i,ind_j) - a(ind_i,ind_j)*k1_til)
        
        fluxes_split(3) = 2.0d0*(gama - 1.0d0)*eigen1*v(ind_i,ind_j) + &
                                eigen3*(v(ind_i,ind_j) + a(ind_i,ind_j)*k2_til) + &
                                eigen4*(v(ind_i,ind_j) - a(ind_i,ind_j)*k2_til)
        
        wII             = (3.0d0 - gama)*(eigen3 + eigen4)*a(ind_i,ind_j)**2.0d0/(2.0d0*(gama-1.0d0))                                   
        
        fluxes_split(4) = (gama - 1.0d0)*eigen1*(u(ind_i,ind_j)**2.0d0+v(ind_i,ind_j)**2.0d0) + &
                          0.50d0*eigen3*((u(ind_i,ind_j)+a(ind_i,ind_j)*k1_til)**2.0d0 + &
                          (v(ind_i,ind_j)+a(ind_i,ind_j)*k2_til)**2.0d0) + &
                          0.50d0*eigen4*((u(ind_i,ind_j)-a(ind_i,ind_j)*k1_til)**2.0d0 + &
                          (v(ind_i,ind_j)-a(ind_i,ind_j)*k2_til)**2.0d0) + wII               

        temp = rho(ind_i,ind_j)*metric_jacobian(ind_i,ind_j)/(2.0d0*gama)
        fluxes_split(1) = temp*fluxes_split(1)
        fluxes_split(2) = temp*fluxes_split(2)
        fluxes_split(3) = temp*fluxes_split(3)
        fluxes_split(4) = temp*fluxes_split(4) 

    end function fluxes_split

    subroutine calculate_fluxes(flux_ksi_pos,flux_ksi_neg,flux_eta_pos,flux_eta_neg)
        implicit none
        real(8),dimension(imax,jmax,dim),intent(out)     :: flux_ksi_pos, flux_ksi_neg
        real(8),dimension(imax,jmax,dim),intent(out)     :: flux_eta_pos, flux_eta_neg
        real(8),dimension(dim)                           :: eig_v_ksi, eig_v_eta
        real(8),dimension(dim)                           :: aux_ksi_pos, aux_ksi_neg
        real(8),dimension(dim)                           :: aux_eta_pos, aux_eta_neg
        real(8),dimension(dim)                           :: eta_neg, eta_pos
        real(8),dimension(dim)                           :: ksi_neg, ksi_pos
        integer(4)                                       :: aux_i, aux_j, pos
        
        flux_eta_neg = 0.0d0
        flux_eta_pos = 0.0d0
        flux_ksi_neg = 0.0d0
        flux_ksi_pos = 0.0d0
    
        do aux_j = 1, jmax
            do aux_i = 1, imax
                
                rho(aux_i,aux_j) = Q_barra(aux_i,aux_j,1)/metric_jacobian(aux_i,aux_j)
                u(aux_i,aux_j) = Q_barra(aux_i,aux_j,2)/Q_barra(aux_i,aux_j,1)
                v(aux_i,aux_j) = Q_barra(aux_i,aux_j,3)/Q_barra(aux_i,aux_j,1)
                p(aux_i,aux_j) = (gama-1.0d0) * (Q_barra(aux_i,aux_j,4)/metric_jacobian(aux_i,aux_j) & 
                        - 0.5d0*rho(aux_i,aux_j)*( (u(aux_i,aux_j))**2.0d0 + (v(aux_i,aux_j))**2.0d0))
                a(aux_i,aux_j) = sqrt(gama*p(aux_i,aux_j)/rho(aux_i,aux_j))           

                if ( aux_j /= 1 .and. aux_j/= jmax .and. aux_i /= 1 .and. aux_i/= imax ) then
                    U_contravariant(aux_i,aux_j) = u(aux_i,aux_j)*ksi_x(aux_i,aux_j) + v(aux_i,aux_j)*ksi_y(aux_i,aux_j)
                    V_contravariant(aux_i,aux_j) = u(aux_i,aux_j)*eta_x(aux_i,aux_j) + v(aux_i,aux_j)*eta_y(aux_i,aux_j)
                endif

                eig_v_ksi = diag_ksi(U_contravariant(aux_i,aux_j),a(aux_i,aux_j),ksi_x(aux_i,aux_j),ksi_y(aux_i,aux_j),dim)
                eig_v_eta = diag_eta(V_contravariant(aux_i,aux_j),a(aux_i,aux_j),eta_x(aux_i,aux_j),eta_y(aux_i,aux_j),dim)
                
                do pos = 1, dim
                    ksi_pos(pos) = 0.5d0*(eig_v_ksi(pos) + abs(eig_v_ksi(pos)))
                    ksi_neg(pos) = 0.5d0*(eig_v_ksi(pos) - abs(eig_v_ksi(pos)))

                    eta_pos(pos) = 0.5d0*(eig_v_eta(pos) + abs(eig_v_eta(pos)))
                    eta_neg(pos) = 0.5d0*(eig_v_eta(pos) - abs(eig_v_eta(pos)))
                end do
                
                aux_ksi_pos = fluxes_split(ksi_x(aux_i,aux_j),ksi_y(aux_i,aux_j),ksi_pos(1),ksi_pos(3),ksi_pos(4),aux_i,aux_j)
                aux_ksi_neg = fluxes_split(ksi_x(aux_i,aux_j),ksi_y(aux_i,aux_j),ksi_neg(1),ksi_neg(3),ksi_neg(4),aux_i,aux_j)
                aux_eta_pos = fluxes_split(eta_x(aux_i,aux_j),eta_y(aux_i,aux_j),eta_pos(1),eta_pos(3),eta_pos(4),aux_i,aux_j)
                aux_eta_neg = fluxes_split(eta_x(aux_i,aux_j),eta_y(aux_i,aux_j),eta_neg(1),eta_neg(3),eta_neg(4),aux_i,aux_j)
                
                flux_ksi_pos(aux_i,aux_j,1:dim) = aux_ksi_pos(1:dim)
                flux_ksi_neg(aux_i,aux_j,1:dim) = aux_ksi_neg(1:dim)
                flux_eta_pos(aux_i,aux_j,1:dim) = aux_eta_pos(1:dim)
                flux_eta_neg(aux_i,aux_j,1:dim) = aux_eta_neg(1:dim)
                
            end do
        end do
        
    end subroutine calculate_fluxes
    
end module fluxes_pos_neg