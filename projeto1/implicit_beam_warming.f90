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
real(8),dimension(:,:),allocatable           :: deltaQ, deltaQ_til
!linear system variables
!real(8),dimension(:,:,:),allocatable           :: Ax_sys, Ay_sys
real(8),dimension(:,:,:),allocatable           :: Bx_sys, By_sys
!
integer index_i, index_j
!
!
!
allocate(A_barra(dim,dim), B_barra(dim,dim), M_barra(dim,dim))
allocate(A_minus(dim,dim,imax), A_plus(dim,dim,imax))
allocate(B_minus(dim,dim,jmax), B_plus(dim,dim,jmax))
allocate(deltaQ_til(dim,imax), deltaQ(dim,jmax))
allocate(Id_x(dim,dim,imax),Id_y(dim,dim,jmax))
allocate(Bx_sys(dim,imax),By_sys(dim,jmax))
!allocate(Ax_sys(imax*dim,imax*dim),Ay_sys(jmax*dim,jmax*dim))
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
!Ax_sys      = 0.0d0
!Ay_sys      = 0.0d0
Bx_sys      = 0.0d0
By_sys      = 0.0d0
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
j_sol = 1
do while ( j_sol <= jmax .and. i_sol <= imax )
    !
    i = 1
    !
            call jacobian_ksi(u(i+1,j_sol),v(i+1,j_sol),Q_barra(i+1,j_sol,4),Q_barra(i+1,j_sol,1),ksi_x(i+1,j_sol),ksi_y(i+1,j_sol),dim,A_barra)
            !
            do index_i = 1, dim
                do index_j = 1, dim
                    A_plus(index_i,index_j,i) = 0.5d0*delta_t(i,j_sol)*A_barra(index_i,index_j)
                end do
            end do
    !
    do i = 2, imax - 1
            !
            call jacobian_ksi(u(i+1,j_sol),v(i+1,j_sol),Q_barra(i+1,j_sol,4),Q_barra(i+1,j_sol,1),ksi_x(i+1,j_sol),ksi_y(i+1,j_sol),dim,A_barra)
            !
            do index_i = 1, dim
                do index_j = 1, dim
                    A_plus(index_i,index_j,i) = 0.5d0*delta_t(i,j_sol)*A_barra(index_i,index_j)
                end do
            end do
            !
            call jacobian_ksi(u(i-1,j_sol),v(i-1,j_sol),Q_barra(i-1,j_sol,4),Q_barra(i-1,j_sol,1),ksi_x(i-1,j_sol),ksi_y(i-1,j_sol),dim,A_barra)
            !
            do index_i = 1, dim
                do index_j = 1, dim
                    A_minus(index_i,index_j,i) = -0.5d0*delta_t(i,j_sol)*A_barra(index_i,index_j)
                end do
            end do
            !
    end do
    !
    i = imax
    !
            call jacobian_ksi(u(i-1,j_sol),v(i-1,j_sol),Q_barra(i-1,j_sol,4),Q_barra(i-1,j_sol,1),ksi_x(i-1,j_sol),ksi_y(i-1,j_sol),dim,A_barra)
            !
            do index_i = 1, dim
                do index_j = 1, dim
                    A_minus(index_i,index_j,i) = -0.5d0*delta_t(i,j_sol)*A_barra(index_i,index_j)
                end do
            end do
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
    !       write (99,209)  (Id_x(1,i,j_sol),j=1,dim)
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
    !
        do i = 1, imax
            call residue(i,j_sol)
            Bx_sys(1,i) = -residue1(i,j_sol)
            Bx_sys(2,i) = -residue2(i,j_sol)
            Bx_sys(3,i) = -residue3(i,j_sol)
            Bx_sys(4,i) = -residue4(i,j_sol)
        end do
    !
    ! we need to solve the block tridiagonal system j times
    ! blktriad(maind,lower,upper,id,md,xb,x)
    call blktriad(Id_x,A_minus,A_plus,dim,imax,Bx_sys,deltaQ_til) 
    !
    ! add one to the index loop
    !
    j_sol = j_sol + 1
    !
    ! *******************
    !
    !
    ! step ii)
    !
    i_sol = 1
    !
    j = 1
    !
    call jacobian_eta(u(i_sol,j+1),v(i_sol,j+1),Q_barra(i_sol,j+1,4),Q_barra(i_sol,j+1,1),eta_x(i_sol,j+1),eta_y(i_sol,j+1),dim,B_barra)
    !
    !
    do index_i = 1, dim
        do index_j = 1, dim
            B_plus(j,index_i,index_j) = 0.5d0*delta_t(i_sol,j)*B_barra(index_i,index_j)
        end do
    end do
    !
    !
    do j = 2, jmax - 1
        !
        !
        call jacobian_eta(u(i_sol,j+1),v(i_sol,j+1),Q_barra(i_sol,j+1,4),Q_barra(i_sol,j+1,1),eta_x(i_sol,j+1),eta_y(i_sol,j+1),dim,B_barra)
        !
        !
        do index_i = 1, dim
            do index_j = 1, dim
                B_plus(j,index_i,index_j) = 0.5d0*delta_t(i_sol,j)*B_barra(index_i,index_j)
            end do
        end do
        !
        !
        call jacobian_eta(u(i_sol,j-1),v(i_sol,j-1),Q_barra(i_sol,j-1,4),Q_barra(i_sol,j-1,1),eta_x(i_sol,j-1),eta_y(i_sol,j-1),dim,B_barra)
        !
        !
        do index_i = 1, dim
            do index_j = 1, dim
                B_minus(j,index_i,index_j) = -0.5d0*delta_t(i_sol,j)*B_barra(index_i,index_j)
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
    call jacobian_eta(u(i_sol,j-1),v(i_sol,j-1),Q_barra(i_sol,j-1,4),Q_barra(i_sol,j-1,1),eta_x(i_sol,j-1),eta_y(i_sol,j-1),dim,B_barra)
    !
    !
    do index_i = 1, dim
        do index_j = 1, dim
            B_minus(j,index_i,index_j) = -0.5d0*delta_t(i_sol,j)*B_barra(index_i,index_j)
        end do
    end do
    !
    ! now create By_sys
    !
    do j = 1, jmax
        By_sys(1,j) = deltaQ_til(1,i_sol)
        By_sys(2,j) = deltaQ_til(2,i_sol)
        By_sys(3,j) = deltaQ_til(3,i_sol)
        By_sys(4,j) = deltaQ_til(4,i_sol)
    end do
    !
    ! solve the block tridiagonal i times
    ! call blktriad(Id_x,A_minus,A_plus,dim,imax,Bx_sys,deltaQ_til) 
    call blktriad(Id_y,B_minus,B_plus,dim,jmax,By_sys,deltaQ)
    !
    ! add one to the index loop
    !
i_sol = i_sol + 1
!
! update Q_barra
!
        Q_barra(i,j,1) = deltaQ(i,j,1) + Q_barra(i,j,1)         
        Q_barra(i,j,2) = deltaQ(i,j,2) + Q_barra(i,j,2)
        Q_barra(i,j,3) = deltaQ(i,j,3) + Q_barra(i,j,3)
        Q_barra(i,j,4) = deltaQ(i,j,4) + Q_barra(i,j,4)   
end do
!
!
!deallocate(Ay_sys)
!deallocate(Ax_sys)
deallocate(Id_x)
deallocate(Bx_sys)
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
