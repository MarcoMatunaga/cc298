subroutine residue_flux_vector_splitting(ind_i,ind_j,Epos,Eneg,Fpos,Fneg,fluxresidue)
    use vars
    implicit none
    integer(4),intent(in)                               :: ind_i, ind_j
    real(8),dimension(imax,jmax,dim),intent(in)         :: Epos, Eneg, Fpos, Fneg
    real(8),dimension(imax,jmax,dim),intent(out)        :: fluxresidue
    !
    max_residue = -100.0d0
    fluxresidue(ind_i,ind_j,1) = delta_t(ind_i,ind_j)*( Epos(ind_i,ind_j,1) - Epos(ind_i-1,ind_j,1) &
                                                       + Eneg(ind_i+1,ind_j,1) - Eneg(ind_i,ind_j,1) &
                                                       + Fpos(ind_i,ind_j,1) - Fpos(ind_i,ind_j-1,1) &
                                                       + Fneg(ind_i,ind_j+1,1) - Fneg(ind_i,ind_j,1) )
    !
    fluxresidue(ind_i,ind_j,2) = delta_t(ind_i,ind_j)*( Epos(ind_i,ind_j,2) - Epos(ind_i-1,ind_j,2) &
                                                       + Eneg(ind_i+1,ind_j,2) - Eneg(ind_i,ind_j,2) &
                                                       + Fpos(ind_i,ind_j,2) - Fpos(ind_i,ind_j-1,2) &
                                                       + Fneg(ind_i,ind_j+1,2) - Fneg(ind_i,ind_j,2) )
    !
    fluxresidue(ind_i,ind_j,3) = delta_t(ind_i,ind_j)*( Epos(ind_i,ind_j,3) - Epos(ind_i-1,ind_j,3) &
                                                       + Eneg(ind_i+1,ind_j,3) - Eneg(ind_i,ind_j,3) &
                                                       + Fpos(ind_i,ind_j,3) - Fpos(ind_i,ind_j-1,3) &
                                                       + Fneg(ind_i,ind_j+1,3) - Fneg(ind_i,ind_j,3) )
    !
    fluxresidue(ind_i,ind_j,4) = delta_t(ind_i,ind_j)*( Epos(ind_i,ind_j,4) - Epos(ind_i-1,ind_j,4) &
                                                       + Eneg(ind_i+1,ind_j,4) - Eneg(ind_i,ind_j,4) &
                                                       + Fpos(ind_i,ind_j,4) - Fpos(ind_i,ind_j-1,4) &
                                                       + Fneg(ind_i,ind_j+1,4) - Fneg(ind_i,ind_j,4) )
    !
    if ( (abs(fluxresidue(ind_i,ind_j,1))) > max_residue ) max_residue = (abs(fluxresidue(ind_i,ind_j,1)))
    if ( (abs(fluxresidue(ind_i,ind_j,2))) > max_residue ) max_residue = (abs(fluxresidue(ind_i,ind_j,2)))
    if ( (abs(fluxresidue(ind_i,ind_j,3))) > max_residue ) max_residue = (abs(fluxresidue(ind_i,ind_j,3)))
    if ( (abs(fluxresidue(ind_i,ind_j,4))) > max_residue ) max_residue = (abs(fluxresidue(ind_i,ind_j,4)))
    !
    max_residue = log10(max_residue)
    if (isnan(max_residue)) stop 'Divergiu :('
    !
    !
end subroutine residue_flux_vector_splitting