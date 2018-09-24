module functions
    use vars
    implicit none

        contains

    real function artificial_dissipation_D4_ksi(art_i,art_j,pos) 
        implicit none
        integer, intent(in) :: art_i, art_j, pos
        double precision :: dis_ksi
        !
        if (art_i == 2 .or. art_i == imax - 1 ) then
        !
                    dis_ksi = -eps_e*( Q_dis(art_i+1,art_j,pos) - 2.0d0*Q_dis(art_i,art_j,pos) &
                                    + Q_dis(art_i-1,art_j,pos) )
                    write(*,*) art_i, dis_ksi
        !
        else
        !
                    dis_ksi = eps_e*( Q_dis(art_i+2,art_j,pos) - 4.0d0*Q_dis(art_i+1,art_j,pos) &
                                    + 6.0d0*Q_dis(art_i,art_j,pos) - 4.0d0*Q_dis(art_i-1,art_j,pos) & 
                                    + Q_dis(art_i-2,art_j,pos) )
        end if
        !
        !
    end function artificial_dissipation_D4_ksi
    !
    !
    !
    real function artificial_dissipation_D4_eta(art_i,art_j,pos) 
        implicit none
        integer, intent(in) :: art_i, art_j, pos
        double precision :: dis_eta
        !
        if (art_j == 2 .or. art_j == jmax - 1 ) then
        !
                    dis_eta = -eps_e*( Q_dis(art_i,art_j+1,pos) - 2.0d0*Q_dis(art_i,art_j,pos) &
                                    + Q_dis(art_i,art_j-1,pos) )
        !
        else
        !
                    dis_eta = eps_e*( Q_dis(art_i,art_j+2,pos) - 4.0d0*Q_dis(art_i,art_j+1,pos) &
                                    + 6.0d0*Q_dis(art_i,art_j,pos) - 4.0d0*Q_dis(art_i,art_j-1,pos) & 
                                    + Q_dis(art_i,art_j-2,pos) )
        end if
        !
        !
    end function artificial_dissipation_D4_eta
    !
    !
    !
    real function artificial_dissipation_D2_eta(art_i,art_j,pos) 
        implicit none
        integer, intent(in) :: art_i, art_j, pos
        double precision :: dis_eta
        !
        ! calculate in everybody 
        !
                    dis_eta = -eps_e*( Q_dis(art_i,art_j+1,pos) - 2.0d0*Q_dis(art_i,art_j,pos) &
                                    + Q_dis(art_i,art_j-1,pos) )
        !
        !
    end function artificial_dissipation_D2_eta
    !
    !
    !
    real function artificial_dissipation_D2_ksi(art_i,art_j,pos) 
        implicit none
        integer, intent(in) :: art_i, art_j, pos
        double precision :: dis_ksi
        !
        ! calculate in everybody 
        !
                    dis_ksi = -eps_e*( Q_dis(art_i+1,art_j,pos) - 2.0d0*Q_dis(art_i,art_j,pos) &
                                    + Q_dis(art_i-1,art_j,pos) )
        !
        !
    end function artificial_dissipation_D2_ksi
    !
    !
    !
end module functions