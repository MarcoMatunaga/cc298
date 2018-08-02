module vars
implicit none
!
!
!
real(8),dimension(:,:,:),allocatable             :: meshx, meshy, meshz
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
allocate(meshx(imax,jmax,kmax), meshy(imax,jmax,kmax))
allocate(meshz(imax,jmax,kmax))
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
subroutine mesh
use vars
implicit none
!
! construct the mesh 
!
do k = 1, kmax
	do j = 1, jmax
		do i = 1, imax
			read(1,*) meshx(i,j,k),meshy(i,j,k),meshz(i,j,k)
		end do
	end do
end do
print *, meshx(3,1,1)
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
			write(2,*) meshx(i,j,k),meshy(i,j,k)
		end do
	end do
end do
close(2)
!
!
!
end subroutine mesh
