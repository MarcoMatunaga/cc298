!
! euler explicit
!
subroutine euler_explicit
use vars
implicit none
!
!
!
do j = 2, jmax - 1
        do i = 2, imax - 1
residue1 = 0.5d0*A*((Q(i+1,j,1) - Q(i-1,j,1))/(meshx(i+1,j) - meshx(i,j))) - 0.5d0*B*((Q(i+1,j,1) -
Q(i-1,j,1))/(meshy(i+1,j) - meshy(i,j)))
Q(i,j,1) = Q(i,j,1) - residue1 
!
residue2 = 0.5d0*A*((Q(i+1,j,2) - Q(i-1,j,2))/(meshx(i+1,j) - meshx(i,j))) - 0.5d0*B*((Q(i+1,j,2) -
Q(i-1,j,2))/(meshy(i+1,j) - meshy(i,j)))
Q(i,j,2) = Q(i,j,2) - residue2  
!
residue3 = 0.5d0*A*((Q(i+1,j,3) - Q(i-1,j,3))/(meshx(i+1,j) - meshx(i,j))) - 0.5d0*B*((Q(i+1,j,3) -
Q(i-1,j,3))/(meshy(i+1,j) - meshy(i,j)))
Q(i,j,3) = Q(i,j,3) - residue3
!
residue4 = 0.5d0*A*((Q(i+1,j,4) - Q(i-1,j,4))/(meshx(i+1,j) - meshx(i,j))) - 0.5d0*B*((Q(i+1,j,4) -
Q(i-1,j,4))/(meshy(i+1,j) - meshy(i,j)))
Q(i,j,4) = Q(i,j,4) - residue4
        end do
end do
!
end subroutine 
