module vars
implicit none
!
!
!
real(8),dimension(:,:),allocatable               :: meshx, meshy	
real(8),dimension(:,:),allocatable  		 :: T, p
real(8)						 :: T_total, p_total
integer(4)				         :: i,j,k
integer(4)					 :: imax,jmax,kmax
!
!
!
end module vars
!
!
!
program proj1
use vars 
implicit none
!
!
!
open(1,file='mesh') 
read(1,*) imax,jmax,kmax
!
!
!
allocate(meshx(imax,jmax), meshy(imax,jmax))
allocate(meshz(imax,jmax))
!
!
call mesh
!
!
deallocate(meshx, meshy, meshz)
end program proj1
!
!
!
subroutine initial_conditions
use vars
implicit none
!
! inlet boundary
! 
! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K
!
T_total = 0.555556d0*531.2d0
P_total = 47.880258888889d0*2117.0d0
!
!
!
u = 0
v = 0
T = T_total
p = p_total

end subroutine initial_conditions
!
!
!
subroutine boundary_conditions
use vars
implicit none
!
! inlet boundary
! 
! 1 lb/ft**2 = 47.880258888889
! 1 R = 0.555556K
!
T_total = 0.555556d0*531.2d0
P_total = 47.880258888889d0*2117.0d0
!
!
!
u = 0
v = 0
T = T_total
p = p_total
end subroutine boundary_conditions
!
!
!
subroutine mesh
use vars
implicit none
!
! construct the mesh 
!
	do j = 1, jmax
		do i = 1, imax
			read(1,*) meshx(i,j),meshy(i,j),meshz(i,j)
		end do
	end do
close(1)
!
! arquivo tecplot
!
open(2,file='mesh.dat')
write(2,*) 'TITLE = "Projeto1" '
write(2,*) 'VARIABLES = "X" "Y" '
WRITE(2,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
do k = 1, kmax
	do j = 1, jmax
		do i = 1, imax
			write(2,*) meshx(i,j),meshy(i,j)
		end do
	end do
end do
close(2)
!
!
!
end subroutine mesh
