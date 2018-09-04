!
! implicit_beam_warming
!
subroutine implicit_beam_warming
    use vars
    implicit none
!
real(8),dimension(:,:),allocatable           :: A_barra, B_barra, M_barra, A_plus, A_minus
real(8),dimension(:,:),allocatable           :: B_plus, B_minus
real(8),dimension(:,:),allocatable           :: Id
real(8),dimension(:),allocatable             :: deltaQ, deltaQ_til, B_sys
real(8),dimension(:,:),allocatable           :: Ax_sys, Ay_sys
!
!
!
allocate(A_barra(dim,dim), B_barra(dim,dim), M_barra(dim,dim))
allocate(A_minus(dim,dim), A_plus(dim,dim))
allocate(B_minus(dim,dim), B_plus(dim,dim))
allocate(deltaQ_til(dim), deltaQ(dim))
allocate(Id(dim,dim))
allocate(B_sys(dim),Ax_sys(imax*dim,imax*dim),Ay_sys(jmax*dim,jmax*dim))
!
! create an identy matrix
!
do j = 1, dim
    do i = 1, dim
        if(i == j) Id(i,j) = 1.0d0 
    end do
end do
!
! initialize the variables for the linear system
!
deltaQ_til  = 0.0d0
deltaQ      = 0.0d0
Ax_sys       = 0.0d0
Ay_sys       = 0.0d0
B_sys       = 0.0d0
!
! residue - pick the maximum residue of them 
! which means the norma_infinity
!
max_residue = -1.0d0
!
! remember that the Q_barra is multiplied by metric_jacobian
! when i call the jacobian i divide energy by rho which 
! cancels the metric_jacobian
!
!***********************************************************
! melhor ter feito uma funcao para as jacobians
!***********************************************************
!
!
! refs pulliam Arc2d and Pulliam fundamental algorithms of
! coputational fluid dynamics
!
! the algorithm consists of two steps 
! in each step you solve a linear system like A X = B
! i) Solve eq 4.119 for deltaQ til
!     
! ii) Solve eq 4.120 for deltaQ barra 
!
! step i) a) I need: A, deltaQ_til, E_Barra, F_barra
!
do j = 2, jmax - 1
    do i = 2, imax - 1
        B_sys(1) = -delta_t(i,j)*0.5d0*( E_Barra(i+1,j,1) - E_Barra(i-1,j,1) + F_barra(i,j+1,1) - F_barra(i,j-1,1) )
        B_sys(2) = -delta_t(i,j)*0.5d0*( E_Barra(i+1,j,2) - E_Barra(i-1,j,2) + F_barra(i,j+1,2) - F_barra(i,j-1,2) )
        B_sys(3) = -delta_t(i,j)*0.5d0*( E_Barra(i+1,j,3) - E_Barra(i-1,j,3) + F_barra(i,j+1,3) - F_barra(i,j-1,3) )
        B_sys(4) = -delta_t(i,j)*0.5d0*( E_Barra(i+1,j,4) - E_Barra(i-1,j,4) + F_barra(i,j+1,4) - F_barra(i,j-1,4) )
    end do
end do
        call jacobian_ksi(u(i+1,j),v(i+1,j),Q_barra(i+1,j,4),Q_barra(i+1,j,1),ksi_x(i+1,j),ksi_y(i+1,j),dim,A_barra)
        A_plus  = A_barra
        call jacobian_ksi(u(i-1,j),v(i-1,j),Q_barra(i-1,j,4),Q_barra(i-1,j,1),ksi_x(i-1,j),ksi_y(i-1,j),dim,A_barra)
        A_minus = A_barra
        call create_block_tridiagonal(A_plus)
!
!
!
call thomas_block_tridiagonal 
!
! step i) b) mount the left matrix of the linear system
!
call jacobian_eta(u(i,j+1),v(i,j+1),Q_barra(i,j+1,4),Q_barra(i,j+1,1),eta_x(i,j+1),eta_y(i,j+1),dim,B_barra)
B_plus  = B_barra
call jacobian_eta(u(i,j-1),v(i,j-1),Q_barra(i,j-1,4),Q_barra(i,j-1,1),eta_x(i,j-1),eta_y(i,j-1),dim,B_barra)
B_minus = B_barra
!
!
!
deallocate(A_barra, B_barra, M_barra)
deallocate(A_minus, A_plus)
deallocate(B_minus, B_plus)
deallocate(deltaQ_til, deltaQ)
deallocate(Id)
deallocate(B_sys,Ax_sys,Ay_sys)
!
!
!
end subroutine implicit_beam_warming
