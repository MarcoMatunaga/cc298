subroutine mesh
use vars
implicit none

! construct the mesh
! the mesh created using the icem software does not have a 
! symmetry line, therefore it is necessary to add one more 
! line and the loop is from one to jmax -1

do j = 1, jmax - 1
    do i = 1, imax
        read(1,*) meshx(i,j),meshy(i,j) 
    end do
end do
!**************************************************
! write the symmetry line CREATE A FUNCTION FOR THIS SITUATION
!**************************************************
j = jmax
do i = 1, imax
    meshx(i,j) = meshx(i,jmax-2)
    meshy(i,j) = meshy(i,jmax-1) + (meshy(i,jmax-1) - meshy(i,jmax-2))
end do

! arquivo tecplot

open(2,file='mesh.dat')
write(2,*) 'TITLE = "Projeto1" '
write(2,*) 'VARIABLES = "X" "Y" '
WRITE(2,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
do j = 1, jmax
    do i = 1, imax
        write(2,*) meshx(i,j), meshy(i,j)
    end do
end do

close(1)

end subroutine mesh
