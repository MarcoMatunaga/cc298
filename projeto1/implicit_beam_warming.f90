!
! implicit_beam_warming
!
subroutine implicit_beam_warming
    use vars
    implicit none
!
real(8),dimension(:,:),allocatable             :: A_barra, B_barra, M_barra
real(8),dimension(:,:,:),allocatable           :: A_plus, A_minus
real(8),dimension(:,:,:),allocatable           :: B_plus, B_minus
real(8),dimension(:,:,:),allocatable           :: Id_x, Id_y
real(8),dimension(:,:,:),allocatable           :: deltaQ, deltaQ_til
! linear system variables
! real(8),dimension(:,:,:),allocatable           :: Ax_sys, Ay_sys
real(8),dimension(:,:,:),allocatable           :: Bx_sys, By_sys
!
integer index_i, index_j
!
!
!
allocate(A_barra(dim,dim), B_barra(dim,dim), M_barra(dim,dim))
allocate(A_minus(imax,dim,dim), A_plus(imax,dim,dim))
allocate(B_minus(jmax,dim,dim), B_plus(jmax,dim,dim))
allocate(deltaQ_til(imax,jmax,dim), deltaQ(imax,jmax,dim))
allocate(Id_x(imax,dim,dim),Id_y(jmax,dim,dim))
allocate(Bx_sys(imax,jmax,dim),By_sys(imax,jmax,dim))
! allocate(Ax_sys(imax*dim,imax*dim),Ay_sys(jmax*dim,jmax*dim))
!
! create the two identies matrix
!
do i = 1, dim
     Id_x(1:imax,i,i) = 1.0d0 
     Id_y(1:jmax,i,i) = 1.0d0 
end do
!
! initialize the variables for the linear system
!
deltaQ_til  = 0.0d0
deltaQ      = 0.0d0
! Ax_sys      = 0.0d0
! Ay_sys      = 0.0d0
Bx_sys      = 0.0d0
By_sys      = 0.0d0
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
! let is start the implicit beam warming
!
! step i) a) I need: A, deltaQ_til, E_Barra, F_barra
!
j = 2
do while ( j <= jmax - 1 )
    !
    i = 1
    !
            call jacobian_ksi(u(i+1,j),v(i+1,j),Q_barra(i+1,j,4),Q_barra(i+1,j,1),ksi_x(i+1,j),ksi_y(i+1,j),dim,A_barra)
            !
            do index_i = 1, dim
                do index_j = 1, dim
                    A_plus(i,index_i,index_j) = 0.5d0*delta_t(i,j)*A_barra(index_i,index_j)
                end do
            end do
    !
    do i = 2, imax - 1
            !
            call jacobian_ksi(u(i+1,j),v(i+1,j),Q_barra(i+1,j,4),Q_barra(i+1,j,1),ksi_x(i+1,j),ksi_y(i+1,j),dim,A_barra)
            !
            do index_i = 1, dim
                do index_j = 1, dim
                    A_plus(i,index_i,index_j) = 0.5d0*delta_t(i,j)*A_barra(index_i,index_j)
                end do
            end do
            !
            call jacobian_ksi(u(i-1,j),v(i-1,j),Q_barra(i-1,j,4),Q_barra(i-1,j,1),ksi_x(i-1,j),ksi_y(i-1,j),dim,A_barra)
            !
            do index_i = 1, dim
                do index_j = 1, dim
                    A_minus(i,index_i,index_j) = -0.5d0*delta_t(i,j)*A_barra(index_i,index_j)
                end do
            end do
            !
    end do
    !
    i = imax
    !
            call jacobian_ksi(u(i-1,j),v(i-1,j),Q_barra(i-1,j,4),Q_barra(i-1,j,1),ksi_x(i-1,j),ksi_y(i-1,j),dim,A_barra)
            !
            do index_i = 1, dim
                do index_j = 1, dim
                    A_minus(i,index_i,index_j) = -0.5d0*delta_t(i,j)*A_barra(index_i,index_j)
                end do
            end do
    !
    ! see the subroutine create_block_tridiagonal for information about the parameters
    !
    ! call create_block_tridiagonal(A_minus,Id_x,A_plus,dim,imax,Ax_sys)
!
! use for test or debug 
!
  !    write (*,200)
  !    open(98,file='./debug/Ax_sys')
  !    do i = 1,dim*imax,dim
  !       write (98,209)  Ax_sys(i:i+dim-1,i:i+dim-1)
  !    end do
  ! 200 format (' Matrix A - block tridiagonal')
  !
  !    open(99,file='./debug/inverse')
  !    do i = 1,dim
  !       write (99,209)  (Id_x(1,i,j),j=1,dim)
  !    end do
  !    close(99)
  !    write (*,201)
  ! 201 format (' Matrix I - block tridiagonal')
    !
    ! now create the vector B_sys, i.e., vector b of the system
    !
    ! *********
    ! i really can dropp the j index because the subroutine call 
    ! is made for each value of the index j
    ! *********
        do i = 2, imax - 1
            Bx_sys(i,j,1) = -delta_t(i,j)*0.5d0*( E_Barra(i+1,j,1) - E_Barra(i-1,j,1) + F_barra(i,j+1,1) - F_barra(i,j-1,1) )
            Bx_sys(i,j,2) = -delta_t(i,j)*0.5d0*( E_Barra(i+1,j,2) - E_Barra(i-1,j,2) + F_barra(i,j+1,2) - F_barra(i,j-1,2) )
            Bx_sys(i,j,3) = -delta_t(i,j)*0.5d0*( E_Barra(i+1,j,3) - E_Barra(i-1,j,3) + F_barra(i,j+1,3) - F_barra(i,j-1,3) )
            Bx_sys(i,j,4) = -delta_t(i,j)*0.5d0*( E_Barra(i+1,j,4) - E_Barra(i-1,j,4) + F_barra(i,j+1,4) - F_barra(i,j-1,4) )
        end do
    !
    ! we need to solve the block tridiagonal system j times
    !
    call thomas_block_tridiagonal(1,j,imax,A_minus,Id_x,A_plus,deltaQ_til,Bx_sys) 
    write(*,*) deltaQ_til
    !
    ! add one to the index loop
    !
    j = j + 1
end do
!
!
deallocate(Id_x)
! deallocate(Ax_sys)
deallocate(Bx_sys)
!
! step ii)
!
i = 2
do while ( i <= imax - 1 )
    !
    j = 1
    !
    call jacobian_eta(u(i,j+1),v(i,j+1),Q_barra(i,j+1,4),Q_barra(i,j+1,1),eta_x(i,j+1),eta_y(i,j+1),dim,B_barra)
    !
    !
    do index_i = 1, dim
        do index_j = 1, dim
            B_plus(j,index_i,index_j) = 0.5d0*delta_t(i,j)
        end do
    end do
    !
    !
    do j = 2, jmax - 1
        !
        !
        call jacobian_eta(u(i,j+1),v(i,j+1),Q_barra(i,j+1,4),Q_barra(i,j+1,1),eta_x(i,j+1),eta_y(i,j+1),dim,B_barra)
        !
        !
        do index_i = 1, dim
            do index_j = 1, dim
                B_plus(j,index_i,index_j) = 0.5d0*delta_t(i,j)
            end do
        end do
        !
        !
        call jacobian_eta(u(i,j-1),v(i,j-1),Q_barra(i,j-1,4),Q_barra(i,j-1,1),eta_x(i,j-1),eta_y(i,j-1),dim,B_barra)
        !
        !
        do index_i = 1, dim
            do index_j = 1, dim
                B_minus(j,index_i,index_j) = -0.5d0*delta_t(i,j)
            end do
        end do
        !
        !
    end do
    !
    !
    j = jmax
    !
    !
    call jacobian_eta(u(i,j-1),v(i,j-1),Q_barra(i,j-1,4),Q_barra(i,j-1,1),eta_x(i,j-1),eta_y(i,j-1),dim,B_barra)
    !
    !
    do index_i = 1, dim
        do index_j = 1, dim
            B_minus(j,index_i,index_j) = -0.5d0*delta_t(i,j)
        end do
    end do
    !
    ! now create By_sys
    !
    do j = 1, jmax
        By_sys(i,j,1) = deltaQ_til(i,j,1)
        By_sys(i,j,2) = deltaQ_til(i,j,2)
        By_sys(i,j,3) = deltaQ_til(i,j,3)
        By_sys(i,j,4) = deltaQ_til(i,j,4)
    end do
    !
    ! solve the block tridiagonal i times
    !
    call thomas_block_tridiagonal(i,1,jmax,B_minus,Id_y,B_plus,deltaQ,By_sys)
    !
    ! add one to the index loop
    !
i = i + 1
end do
!
!
! deallocate(Ay_sys)
deallocate(By_sys)
!
!
!
deallocate(A_barra, B_barra, M_barra)
deallocate(A_minus, A_plus)
deallocate(B_minus, B_plus)
deallocate(deltaQ_til, deltaQ)
deallocate(Id_y)
!
!
!
209 format (500f12.6)
end subroutine implicit_beam_warming
