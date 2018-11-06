module fluxes_pos_neg
    use vars
    use diagonalization
    !
    contains
    !
    subroutine eigen_values_calculate(eigen,eigen_pos,eigen_neg)
        implicit none
        real(8),intent(in)                 :: eigen
        real(8),intent(out)                :: eigen_pos,eigen_neg
        !
        !
                eigen_pos = 0.50d0*(eigen + abs(eigen))
                eigen_neg = 0.50d0*(eigen - abs(eigen))
        !
        !
    end subroutine eigen_values_calculate
    !
    function fluxes_split(k1,k2,eigen1,eigen3,eigen4,ind_i,ind_j)
        implicit none
        real(8),intent(in)                  :: k1,k2
        real(8),intent(in)                  :: eigen1,eigen3,eigen4
        integer(4),intent(in)               :: ind_i,ind_j
        real(8),dimension(dim)              :: fluxes_split
        real(8)                             :: rho_temp, u_temp, v_temp, temp
        real(8)                             :: k1_til, k2_til, wII
        !
        !
        rho_temp = Q_barra(ind_i,ind_j,1)/metric_jacobian(ind_i,ind_j)
        u_temp = Q_barra(ind_i,ind_j,2)/Q_barra(ind_i,ind_j,1)
        v_temp = Q_barra(ind_i,ind_j,3)/Q_barra(ind_i,ind_j,1)
        k1_til = k1/sqrt(k1**2.0d0 + k2**2.0d0)
        k2_til = k2/sqrt(k1**2.0d0 + k2**2.0d0)
        !
        fluxes_split(1) = 2.0d0*(gama - 1.0d0)*eigen1 + eigen3 + eigen4
        fluxes_split(2) = 2.0d0*(gama - 1.0d0)*eigen1*u_temp + &
                                   eigen3*(u_temp + a(ind_i,ind_j)*k1_til) + &
                                   eigen4*(u_temp - a(ind_i,ind_j)*k1_til)
        fluxes_split(3) = 2.0d0*(gama - 1.0d0)*eigen1*v_temp + &
                                   eigen3*(v_temp + a(ind_i,ind_j)*k1_til) + &
                                   eigen4*(v_temp - a(ind_i,ind_j)*k1_til)
        wII                      = (3.0d0 - gama)*(eigen3 + eigen4)*a(ind_i,ind_j)**2.0d0/(2.0d0*(gama-1.0d0))                                   
        fluxes_split(4) = (gama - 1.0d0)*eigen1*(u_temp**2.0d0+v_temp**2.0d0) + &
                                   0.50d0*eigen3*((u_temp+a(ind_i,ind_j)*k1_til)**2.0d0 + &
                                   (v_temp+a(ind_i,ind_j)*k2_til)**2.0d0) + &
                                   0.50d0*eigen4*((u_temp-a(ind_i,ind_j)*k1_til)**2.0d0 + &
                                   (v_temp-a(ind_i,ind_j)*k2_til)**2.0d0) + wII               
        !
        temp = rho_temp*metric_jacobian(ind_i,ind_j)/(2.0d0*gama)
        temp = rho_temp*metric_jacobian(ind_i,ind_j)/(2.0d0*gama)
        temp = rho_temp*metric_jacobian(ind_i,ind_j)/(2.0d0*gama)
        temp = rho_temp*metric_jacobian(ind_i,ind_j)/(2.0d0*gama)
        fluxes_split(1) = temp*fluxes_split(1)
        fluxes_split(2) = temp*fluxes_split(2)
        fluxes_split(3) = temp*fluxes_split(3)
        fluxes_split(4) = temp*fluxes_split(4) 
        !
        !
    end function fluxes_split
    !
    subroutine calculate_fluxes(flux_ksi_pos,flux_ksi_neg,flux_eta_pos,flux_eta_neg)
        implicit none
        real(8),dimension(imax,jmax,dim),intent(out)     :: flux_ksi_pos, flux_ksi_neg
        real(8),dimension(imax,jmax,dim),intent(out)     :: flux_eta_pos, flux_eta_neg
        real(8),dimension(dim)                           :: eig_v_ksi, eig_v_eta
        real(8),dimension(dim)                           :: aux_pos_ksi, aux_neg_ksi
        real(8),dimension(dim)                           :: aux_pos_eta, aux_neg_eta
        real(8),dimension(dim)                           :: pos_eta, pos_ksi
        real(8),dimension(dim)                           :: neg_eta, neg_ksi
        integer(4)                                       :: aux_i, aux_j, pos
        !
        do aux_j = 1, jmax
            do aux_i = 1, imax
                !
                eig_v_ksi = diag_ksi(U_contravariant(aux_i,aux_j),a(aux_i,aux_j),ksi_x(aux_i,aux_j),ksi_y(aux_i,aux_j),dim)
                eig_v_eta = diag_eta(V_contravariant(aux_i,aux_j),a(aux_i,aux_j),eta_x(aux_i,aux_j),eta_y(aux_i,aux_j),dim)
                do pos = 1, dim
                    call eigen_values_calculate(eig_v_ksi(pos),pos_ksi(pos),neg_ksi(pos))
                    call eigen_values_calculate(eig_v_eta(pos),pos_eta(pos),neg_eta(pos))
                end do
                aux_pos_ksi = fluxes_split(ksi_x(aux_i,aux_j),ksi_y(aux_i,aux_j),pos_ksi(1),pos_ksi(3),pos_ksi(4),aux_i,aux_j)
                aux_neg_ksi = fluxes_split(ksi_x(aux_i,aux_j),ksi_y(aux_i,aux_j),neg_ksi(1),neg_ksi(3),neg_ksi(4),aux_i,aux_j)
                aux_pos_eta = fluxes_split(eta_x(aux_i,aux_j),eta_y(aux_i,aux_j),pos_eta(1),pos_eta(3),pos_eta(4),aux_i,aux_j)
                aux_neg_eta = fluxes_split(eta_x(aux_i,aux_j),eta_y(aux_i,aux_j),neg_eta(1),neg_eta(3),neg_eta(4),aux_i,aux_j)
                !
                flux_ksi_pos(aux_i,aux_j,1:dim) = aux_pos_ksi(1:dim)
                flux_ksi_neg(aux_i,aux_j,1:dim) = aux_neg_ksi(1:dim)
                flux_eta_pos(aux_i,aux_j,1:dim) = aux_pos_eta(1:dim)
                flux_eta_neg(aux_i,aux_j,1:dim) = aux_neg_eta(1:dim)
                !
            end do
        end do
        !
    end subroutine calculate_fluxes
    !
end module fluxes_pos_neg