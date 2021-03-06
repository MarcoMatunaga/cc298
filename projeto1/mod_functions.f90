module functions
    use vars
    implicit none
    !
        contains
    !
    real(8) function dis_ksi4(art_i,art_j,pos,eps) 
        implicit none
        integer(4), intent(in) :: art_i, art_j, pos
        real(8),    intent(in) :: eps
        !
        if (art_i == 2 .or. art_i == imax - 1) then
        !
                    dis_ksi4 = -eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*( Q_dis(art_i+1,art_j,pos) &
                               - 2.0d0*Q_dis(art_i,art_j,pos) &
                               + Q_dis(art_i-1,art_j,pos) )
        !
        else
        !
                    dis_ksi4 = eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*( Q_dis(art_i+2,art_j,pos) &
                               - 4.0d0*Q_dis(art_i+1,art_j,pos) &
                               + 6.0d0*Q_dis(art_i,art_j,pos) - 4.0d0*Q_dis(art_i-1,art_j,pos) & 
                               + Q_dis(art_i-2,art_j,pos) )
        end if
        !
        !
    end function dis_ksi4
    !
    !
    real(8) function dis_eta4(art_i,art_j,pos,eps) 
        implicit none
        integer(4), intent(in) :: art_i, art_j, pos
        real(8),    intent(in) :: eps
        !
        if (art_j == 2 .or. art_j == jmax - 1) then
        !
                    dis_eta4 = -eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*( Q_dis(art_i,art_j+1,pos) &
                               - 2.0d0*Q_dis(art_i,art_j,pos) &
                               + Q_dis(art_i,art_j-1,pos) )
        !
        else
        !
                    dis_eta4 = eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*( Q_dis(art_i,art_j+2,pos) &
                               - 4.0d0*Q_dis(art_i,art_j+1,pos) &
                               + 6.0d0*Q_dis(art_i,art_j,pos) - 4.0d0*Q_dis(art_i,art_j-1,pos) & 
                               + Q_dis(art_i,art_j-2,pos) )
        end if
        !
        !
    end function dis_eta4
    !
    !
    real(8) function dis_eta2(art_i,art_j,pos,eps) 
        implicit none
        integer(4), intent(in) :: art_i, art_j, pos
        real(8),    intent(in) :: eps
        !
        !
                    dis_eta2 = -eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*( Q_dis(art_i,art_j+1,pos) &
                               - 2.0d0*Q_dis(art_i,art_j,pos) + Q_dis(art_i,art_j-1,pos) )
        !
        !
    end function dis_eta2
    !
    !
    real(8) function dis_ksi2(art_i,art_j,pos,eps) 
        implicit none
        integer(4), intent(in) :: art_i, art_j, pos
        real(8),    intent(in) :: eps
        !
        !
                    dis_ksi2 = -eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*( Q_dis(art_i+1,art_j,pos) &
                               - 2.0d0*Q_dis(art_i,art_j,pos) + Q_dis(art_i-1,art_j,pos) )
        !
        !
    end function dis_ksi2
    !
    !
    real(8) function dis_imp_ksi(art_i,art_j,eps,pos)
        implicit none
        integer(4), intent(in) :: art_i, art_j, pos
        real(8),    intent(in) :: eps
        real(8),dimension(3)   :: vector
        !
        !
                if (pos == 1) then
                    vector(1) = -eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*(1.0d0/metric_jacobian(art_i-1,art_j))
                else if (pos == 2) then 
                    vector(2) = eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*(2.0d0/metric_jacobian(art_i,art_j))
                else
                    vector(3) = -eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*(1.0d0/metric_jacobian(art_i+1,art_j))      
                end if
                    !
                    dis_imp_ksi = vector(pos)
                    !              
        !
        !
    end function dis_imp_ksi
    !
    !
    real(8) function dis_imp_eta(art_i,art_j,eps,pos)
        implicit none
        integer(4), intent(in) :: art_i, art_j, pos
        real(8),    intent(in) :: eps
        real(8),dimension(3)   :: vector
        !
        !
                if (pos == 1) then
                    vector(1) = -eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*(1.0d0/metric_jacobian(art_i,art_j-1))
                else if (pos == 2) then 
                    vector(2) = eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*(2.0d0/metric_jacobian(art_i,art_j))
                else
                    vector(3) = -eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)*(1.0d0/metric_jacobian(art_i,art_j+1)) 
                end if
                    !
                    dis_imp_eta = vector(pos)
                    !
        !
        !
    end function dis_imp_eta
    !
    !
    real(8) function  non_lin_dis_ksi(art_i,art_j,eps2,eps4,pos)
        implicit none
        integer(4), intent(in) :: art_i, art_j, pos
        real(8),    intent(in) :: eps2, eps4
        !
        !
            non_lin_dis_ksi = sigma_ksi(art_i+1,art_j)*metric_jacobian(art_i+1,art_j) &
                              + sigma_ksi(art_i,art_j)*metric_jacobian(art_i,art_j)
            non_lin_dis_ksi = non_lin_dis_ksi*( eps2*(Q_dis(art_i+1,art_j,pos) - Q_dis(art_i,art_j,pos)) &
                              - eps4*(Q_dis(art_i+2,art_j,pos) - 3.0d0*Q_dis(art_i+1,art_j,pos) &
                              + 3.0d0*Q_dis(art_i,art_j,pos) - Q_dis(art_i-1,art_j,pos)))
        !
        !
    end function non_lin_dis_ksi
    !
    !
    reaL(8) function non_lin_dis_eta(art_i,art_j,eps2,eps4,pos)
        implicit none
        integer(4), intent(in) :: art_i, art_j, pos
        real(8),    intent(in) :: eps2, eps4
        !
        !
            non_lin_dis_eta = sigma_eta(art_i,art_j+1)*metric_jacobian(art_i,art_j+1) &
                              + sigma_eta(art_i,art_j)*metric_jacobian(art_i,art_j)
            non_lin_dis_eta = non_lin_dis_eta*( eps2*(Q_dis(art_i,art_j+1,pos) - Q_dis(art_i,art_j,pos)) &
                              - eps4*(Q_dis(art_i,art_j+2,pos) - 3.0d0*Q_dis(art_i,art_j+1,pos) &
                              + 3.0d0*Q_dis(art_i,art_j,pos) - Q_dis(art_i,art_j-1,pos)))
        !
        !
    end function non_lin_dis_eta
    !
    !
    function dis_ksi4_imp(art_i,art_j,pos,eps) 
        implicit none
        integer(4), intent(in) :: art_i, art_j, pos
        real(8),    intent(in) :: eps
        real(8),dimension(5)   :: dis_ksi4_imp
        real(8)                :: mult_factor
        
        if (art_i == 2 .or. art_i == imax - 1) then
        
                    mult_factor = -eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)

                    dis_ksi4_imp(1) = 0.0d0
                    dis_ksi4_imp(2) = mult_factor*(1.0d0/metric_jacobian(art_i+1,art_j))*Q_dis(art_i+1,art_j,pos)
                    dis_ksi4_imp(3) = mult_factor*(1.0d0/metric_jacobian(art_i,art_j))*(-2.0d0)*Q_dis(art_i,art_j,pos)
                    dis_ksi4_imp(4) = mult_factor*(1.0d0/metric_jacobian(art_i-1,art_j))*Q_dis(art_i-1,art_j,pos)
                    dis_ksi4_imp(5) = 0.0d0

        else
        
                    mult_factor = eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)

                    dis_ksi4_imp(1) = mult_factor*(1.0d0/metric_jacobian(art_i+2,art_j))*Q_dis(art_i+2,art_j,pos)
                    dis_ksi4_imp(2) = mult_factor*(1.0d0/metric_jacobian(art_i+1,art_j))*(-4.0d0)*Q_dis(art_i+1,art_j,pos)
                    dis_ksi4_imp(3) = mult_factor*(1.0d0/metric_jacobian(art_i,art_j))*6.0d0*Q_dis(art_i,art_j,pos)
                    dis_ksi4_imp(4) = mult_factor*(1.0d0/metric_jacobian(art_i-1,art_j))*(-4.0d0)*Q_dis(art_i-1,art_j,pos)
                    dis_ksi4_imp(5) = mult_factor*(1.0d0/metric_jacobian(art_i-2,art_j))*Q_dis(art_i-2,art_j,pos)

        end if

    end function dis_ksi4_imp

    function dis_eta4_imp(art_i,art_j,pos,eps) 
        implicit none
        integer(4), intent(in) :: art_i, art_j, pos
        real(8),    intent(in) :: eps
        real(8), dimension(5)  :: dis_eta4_imp
        real(8)                :: mult_factor
        
        if (art_j == 2 .or. art_j == jmax - 1) then
        
                    mult_factor = -eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)

                    dis_eta4_imp(1) = 0.0d0
                    dis_eta4_imp(2) = mult_factor*(1.0d0/metric_jacobian(art_i,art_j+1))*Q_dis(art_i,art_j+1,pos)
                    dis_eta4_imp(3) = mult_factor*(1.0d0/metric_jacobian(art_i,art_j))*(-2.0d0)*Q_dis(art_i,art_j,pos)
                    dis_eta4_imp(4) = mult_factor*(1.0d0/metric_jacobian(art_i,art_j-1))*Q_dis(art_i,art_j-1,pos) 
                    dis_eta4_imp(5) = 0.0d0
        
        else
        
                    mult_factor = eps*delta_t(art_i,art_j)*metric_jacobian(art_i,art_j)

                    dis_eta4_imp(1) = mult_factor*(1.0d0/metric_jacobian(art_i,art_j+2))*Q_dis(art_i,art_j+2,pos)
                    dis_eta4_imp(2) = mult_factor*(1.0d0/metric_jacobian(art_i,art_j+1))*(-4.0d0)*Q_dis(art_i,art_j+1,pos)
                    dis_eta4_imp(3) = mult_factor*(1.0d0/metric_jacobian(art_i,art_j))*6.0d0*Q_dis(art_i,art_j,pos)
                    dis_eta4_imp(4) = mult_factor*(1.0d0/metric_jacobian(art_i,art_j-1))*(-4.0d0)*Q_dis(art_i,art_j-1,pos)
                    dis_eta4_imp(5) = mult_factor*(1.0d0/metric_jacobian(art_i,art_j-2))*Q_dis(art_i,art_j-2,pos)

        end if
        !
        !
    end function dis_eta4_imp

end module functions