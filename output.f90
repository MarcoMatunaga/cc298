subroutine output
use vars 
implicit none
real(8),dimension(:,:),allocatable           :: p_out, u_out, v_out, rho
!
!
allocate(p_out(imax,jmax), u_out(imax,jmax), v_out(imax,jmax), rho(imax,jmax) )
!
! testar as metricas de transformacao
!
! open(3,file='teste.dat')
! write(3,*) 'TITLE = "Projeto1" '
! write(3,*) 'VARIABLES = "X" "Y" "X_ksi" "X_eta" "Y_ksi" "Y_eta" "J" '
! write(3,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
! do j = 1, jmax
!     do i = 1, imax
!         !write(3,'(7es11.3e2)') meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
!         !write(3,*) meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
!         write(3,'(7EN20.10)') meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
!     end do
! end do
! close(3)
!
 do j = 1, jmax
     do i = 1, imax
         rho(i,j)   = Q_barra(i,j,1)/metric_jacobian(i,j)
         u_out(i,j) = Q_barra(i,j,2)/Q_barra(i,j,1)
         v_out(i,j) = Q_barra(i,j,3)/Q_barra(i,j,1)
         p_out(i,j) = (gama - 1.0d0)*((Q_barra(i,j,4)/metric_jacobian(i,j)) &
                      - 0.50d0*( (u_out(i,j)/rho(i,j))**2.0d0 + (v_out(i,j)/rho(i,j))**2.0d0 ))
         !print *, i,j, p_out(i,j), u_out(i,j), v_out(i,j)
     end do
 end do
!
! condicoes de contorno e iniciais
!
open(3,file='teste_boundary.dat')
write(3,*) 'TITLE = "Projeto1" '
write(3,*) 'VARIABLES = "X" "Y" "u" "v" "rho" "p" '
write(3,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
do j = 1, jmax
    do i = 1, imax
        !write(3,'(7es11.3e2)') meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
        !write(3,*) meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
        write(3,'(6EN20.10)') meshx(i,j), meshy(i,j), u_out(i,j), v_out(i,j), rho(i,j), p_out(i,j) 
    end do
end do
close(3)
!
!
!
deallocate(u_out, v_out, p_out, rho)
!
!
!
end subroutine output