subroutine residue_flux_vector_splitting_2nd(ii,jj,Epos,Eneg,Fpos,Fneg,fluxresidue)
    use vars
    implicit none
    integer(4),intent(in)                                   :: ii,jj
    real(8),dimension(imax,jmax,dim),intent(in)             :: Epos, Eneg, Fpos, Fneg
    real(8),dimension(imax,jmax,dim),intent(out)            :: fluxresidue
    integer(4)                                              :: ind_p
    
    max_residue = -100.0d0
    !*******
    !create a function an operator backward
    !*******
    do ind_p = 1, dim
      fluxresidue(ii,jj,ind_p) = delta_t(ii,jj)*0.50d0*( 3.0d0*Epos(ii,jj,ind_p) - 4.0d0*Epos(ii-1,jj,ind_p) &
                                                        + Epos(ii-2,jj,ind_p) &
                                                        + 3.0d0*Fpos(ii,jj,ind_p) - 4.0d0*Fpos(ii,jj-1,ind_p) &
                                                        + Fpos(ii,jj-2,ind_p) & 
                                                        -3.0d0*Eneg(ii,jj,ind_p) + 4.0d0*Eneg(ii+1,jj,ind_p) &
                                                        -Eneg(ii+2,jj,ind_p) &
                                                        -3.0d0*Fneg(ii,jj,ind_p) + 4.0d0*Fneg(ii,jj+1,ind_p) &
                                                        -Fneg(ii,jj+2,ind_p) )
    end do
    
    if ( (abs(fluxresidue(ii,jj,1))) > max_residue ) max_residue = (abs(fluxresidue(ii,jj,1)))
    if ( (abs(fluxresidue(ii,jj,2))) > max_residue ) max_residue = (abs(fluxresidue(ii,jj,2)))
    if ( (abs(fluxresidue(ii,jj,3))) > max_residue ) max_residue = (abs(fluxresidue(ii,jj,3)))
    if ( (abs(fluxresidue(ii,jj,4))) > max_residue ) max_residue = (abs(fluxresidue(ii,jj,4)))
    
    max_residue = log10(max_residue)
    if (isnan(max_residue)) stop 'Divergiu :('

end subroutine residue_flux_vector_splitting_2nd