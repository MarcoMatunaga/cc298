subroutine flux_curvilinear_transformation_2D(V1,x_index_max,y_index_max,dim,eta_x,ksi_x,eta_y,ksi_y,jacobian,V1_barra)
implicit none
real(8), dimension(x_index_max,y_index_max,dim), intent(in)             :: V1                            ! input
real(8), dimension(x_index_max,y_index_max), intent(in)                 :: jacobian                      ! input
integer(4), intent(in)                                                  :: x_index_max,y_index_max,dim   ! input
real(8), dimension(x_index_max,y_index_max,dim),intent(out)             :: V1_barra                      ! output
integer(4)                                                              :: i,j
!
!
! read book 
! Euler equations 
!   - vectors of conserved variables is denoted by Q
!   - fluxes vector are denoted by E and F, in the x and y direction
!     respectively
! V1 = Q, E, F
!
!
! pressure gradient in the curvilinear space
!

!
! temperature gradient in the curvilinear space
!

!
! velocities components in the curvilinear space
!
end subroutine flux_curvilinear_transformation_2D