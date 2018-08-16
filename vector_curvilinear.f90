subroutine curvilinear_transformation_2D(V1,V2,V3,x_index_max,y_index_max,dim,V1_barra,V2_barra,V3_barra)
implicit none
real(8), dimension(:,:,:), intent(in)         :: V1, V2, V3                   ! input
real(8), dimension(:,:,:), intent(out)        :: V1_barra,V2_barra,V3_barra   ! output
integer                                       :: i,j
!
! Euler equations 
!   - vectors of conserved variables is denoted by Q
!   - fluxes vector are denoted by E and F, in the x and y direction
!     respectively
! V1 = Q
! V2 = E
! V3 = F
!
end subroutine curvilinear_transformation_2D