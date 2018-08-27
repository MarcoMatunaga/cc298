program serie2
    ! pontos da malha - usar no maximo 50 pontos
    integer(4)             :: j, jmiddle, jmax
    integer(4)             :: i
    real(8),dimension(50)  :: u, f, u_init
!
! declare the values
!
u = 0.0d0
jmax = 10
jmiddle = jmax/2
print *, jmiddle
!
!
!
!
! initial conditions
!
do i = 1, jmiddle - 1
    u_init(i) = 1.0d0
    u(i) = u_init(i)
end do
u_init(jmiddle) = 0.0d0
u(jmiddle) = u_init(jmiddle)
do i = jmiddle + 1, jmax
    u_init(i) = -1.0d0
    u(i) = u_init(i)
end do
!
! artificial dissipation - second differences
!
do i = 1, jmax
f(i) = - 0.25d0*( u(i+2) - 4.0d0*u(i+1) + 6.0d0*u(i) - 4.0d0*u(i-1) + u(j-2) )
end do
!
! time marching
!
do i = 2, jmax - 1
u(i) = u(i) - 0.50d0*(u(i+1) - u(i-1)) + f(i)
end do
!
!
!
open(1,file="diferenças_segundas")
do i = 1, jmax
write(1,*) i, u_init(i), u(i)
end do
close(1)
!*****************************************************************!
!*****************************************************************!
!*****************************************************************!
!*****************************************************************!
!
! initial conditions
!
do i = 1, jmiddle - 1
    u_init(i) = 1.0d0
    u(i) = u_init(i)
end do
u_init(jmiddle) = 0.0d0
u(jmiddle) = u_init(jmiddle)
do i = jmiddle + 1, jmax
    u_init(i) = -1.0d0
    u(i) = u_init(i)
end do
!
! artificial dissipation - fourth differences
!
do i = 1, jmax
    f(i) = 0.50d0*( u(i+1) - 2.0d0*u(i) + u(i-1) )
end do
!
! time marching
!
do i = 2, jmax - 1
u(i) = u(i) - 0.50d0*(u(i+1) - u(i-1)) + f(i)
end do
!
! output
!
open(2,file="diferenças_quartas")
do i = 1, jmax
write(2,*) i, u_init(i),u(i)
end do
close(2)
end program serie2