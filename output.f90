subroutine output
use vars 
implicit none
!
!
!
open(3,file='teste.dat')
write(3,*) 'TITLE = "Projeto1" '
write(3,*) 'VARIABLES = "X_ksi" "X_eta" "Y_ksi" "Y_eta" "J" '
write(3,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
do j = 1, jmax
    do i = 1, imax
        write(3,'(7es11.3e2)') x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
    end do
end do
close(3)
!
!
!
end subroutine output
