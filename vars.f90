module vars
implicit none
! variaveis do fluido e do escoamento	
real(8),dimension(:,:),allocatable  		 :: T, p, u, v, rho, q_vel, a
real(8),dimension(:,:,:),allocatable 		 :: Q 
real(8)						 :: T_total, p_total
! constantes do fluido
real(8)						 :: gama, c_v, R
! constantes do escoamento 
real(8)						 :: theta 
! constantes matematicas
real(8)						 :: pi, dummy
! indices dos vetores --> geometria, malha
real(8),dimension(:,:),allocatable               :: meshx, meshy
integer(4)				         :: i,j,k
integer(4)					 :: imax,jmax,kmax
! variaveis numericas 
!
!
!
end module vars
