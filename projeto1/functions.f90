module functions
    use vars
    implicit none

        contains

    real function dis_ksi4(art_i,art_j,pos) 
        implicit none
        integer, intent(in) :: art_i, art_j, pos
        !
        write(*,*) Q_dis(art_i+1,art_j,pos), Q_dis(art_i,art_j,pos),  Q_dis(art_i-1,art_j,pos)
        if (art_i == 2 .or. art_i == imax - 1 .or. art_j == 2 .or. art_j == jmax - 1) then
        !
                    dis_ksi4 = -eps_e*( Q_dis(art_i+1,art_j,pos) - 2.0d0*Q_dis(art_i,art_j,pos) &
                               + Q_dis(art_i-1,art_j,pos) )
        !
        else
        !
                    dis_ksi4 = eps_e*( Q_dis(art_i+2,art_j,pos) - 4.0d0*Q_dis(art_i+1,art_j,pos) &
                               + 6.0d0*Q_dis(art_i,art_j,pos) - 4.0d0*Q_dis(art_i-1,art_j,pos) & 
                               + Q_dis(art_i-2,art_j,pos) )
        end if
        !
        !
    end function dis_ksi4
    !
    !
    !
    real function dis_eta4(art_i,art_j,pos) 
        implicit none
        integer, intent(in) :: art_i, art_j, pos
        !
        if (art_i == 2 .or. art_i == imax - 1 .or. art_j == 2 .or. art_j == jmax - 1) then
        !
                    dis_eta4 = -eps_e*( Q_dis(art_i,art_j+1,pos) - 2.0d0*Q_dis(art_i,art_j,pos) &
                               + Q_dis(art_i,art_j-1,pos) )
        !
        else
        !
                    dis_eta4 = eps_e*( Q_dis(art_i,art_j+2,pos) - 4.0d0*Q_dis(art_i,art_j+1,pos) &
                               + 6.0d0*Q_dis(art_i,art_j,pos) - 4.0d0*Q_dis(art_i,art_j-1,pos) & 
                               + Q_dis(art_i,art_j-2,pos) )
        end if
        !
        !
    end function dis_eta4
    !
    !
    !
    real function dis_eta2(art_i,art_j,pos) 
        implicit none
        integer, intent(in) :: art_i, art_j, pos
        !
        !
                    dis_eta2 = -eps_e*( Q_dis(art_i,art_j+1,pos) - 2.0d0*Q_dis(art_i,art_j,pos) &
                                    + Q_dis(art_i,art_j-1,pos) )
        !
        !
    end function dis_eta2
    !
    !
    !
    real function dis_ksi2(art_i,art_j,pos) 
        implicit none
        integer, intent(in) :: art_i, art_j, pos
        !
        !
                    dis_ksi2 = -eps_e*( Q_dis(art_i+1,art_j,pos) - 2.0d0*Q_dis(art_i,art_j,pos) &
                                    + Q_dis(art_i-1,art_j,pos) )
        !
        !
    end function dis_ksi2
    !
    !
    !
end module functions