
! implicit_beam_warming

subroutine implicit_beam_warming
    use vars
    use functions
    use diagonalization
    use functionsDerivatives
    implicit none

real(8),dimension(:,:),allocatable             :: A_barra, B_barra, M_barra
real(8),dimension(:,:,:),allocatable           :: A_plus, A_minus
real(8),dimension(:,:,:),allocatable           :: B_plus, B_minus
real(8),dimension(:,:,:),allocatable           :: Id_x, Id_y
real(8),dimension(:,:,:),allocatable           :: main_x, main_y
real(8),dimension(:,:),allocatable             :: deltaQ_til_i, deltaQ_til_j
real(8),dimension(:,:,:),allocatable           :: deltaQ, deltaQ_til
!linear system variables
!real(8),dimension(:,:,:),allocatable           :: Ax_sys, Ay_sys
real(8),dimension(:,:),allocatable             :: Bx_sys, By_sys
real(8),dimension(:,:),allocatable             :: Identy
real(8)                                        :: L_ksi, L_eta
real(8)                                        :: u, v

integer index_i, index_j, i_sol, j_sol

allocate(A_barra(dim,dim), B_barra(dim,dim), M_barra(dim,dim))
allocate(A_minus(dim,dim,imax), A_plus(dim,dim,imax))
allocate(B_minus(dim,dim,jmax), B_plus(dim,dim,jmax))
allocate(deltaQ_til_i(dim,imax), deltaQ_til_j(dim,jmax))
allocate(deltaQ_til(dim,imax,jmax), deltaQ(imax,jmax,dim))
allocate(Id_x(dim,dim,imax),Id_y(dim,dim,jmax))
allocate(main_x(dim,dim,imax),main_y(dim,dim,jmax))
allocate(Bx_sys(dim,imax),By_sys(dim,jmax))
allocate(Identy(dim,dim))

! debug files
open(10, file='debug_fdJAC') 
open(11, file='JAC_correct') 

! create the two identies matrix

!************
! we can eliminate the variables in the y direction
! it is possible to allocate and deallocate the arrays
!************
Id_x = 0.0d0
Id_y = 0.0d0
do i = 1, dim
     Id_x(i,i,1:imax) = 1.0d0 
     Id_y(i,i,1:jmax) = 1.0d0
     Identy(i,i) = 1.0d0 
end do
!
! initialize the variables for the linear system
!
deltaQ_til_i  = 0.0d0
deltaQ_til_j  = 0.0d0
deltaQ      = 0.0d0
main_x = 0.0d0
main_y = 0.0d0
Bx_sys      = 0.0d0
By_sys      = 0.0d0
!
!
do j = 1, jmax
        do i = 1, imax
            Q_dis(i,j,1) = Q_barra(i,j,1)/metric_jacobian(i,j)
            Q_dis(i,j,2) = Q_barra(i,j,2)/metric_jacobian(i,j)
            Q_dis(i,j,3) = Q_barra(i,j,3)/metric_jacobian(i,j)
            Q_dis(i,j,4) = Q_barra(i,j,4)/metric_jacobian(i,j)
        end do
end do
!
! remember that the Q_barra is multiplied by metric_jacobian
! when i call the jacobian i divide energy by rho which 
! cancels the metric_jacobian
!
! ***********************************************************
!   melhor ter feito uma funcao para as jacobians
! ***********************************************************
!
! refs pulliam Arc2d and Pulliam fundamental algorithms of
! computational fluid dynamics
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
i_sol = 2 ! change j_sol by j
j_sol = 2 ! idem for i_sol
!looping_j_sol: do while ( j_sol <= jmax - 1 )
looping_j_sol: do j_sol = 2, jmax - 1
    
    do i = 2, imax - 1
            
            L_ksi = dis_imp_ksi(i,j_sol,eps_dis_i,3)
            u = Q_barra(i+1,j_sol,2)/Q_barra(i+1,j_sol,1)
            v = Q_barra(i+1,j_sol,3)/Q_barra(i+1,j_sol,1)
call jacobian_ksi(u,v,Q_barra(i+1,j_sol,4),Q_barra(i+1,j_sol,1),ksi_x(i+1,j_sol),ksi_y(i+1,j_sol),dim,A_barra)
call FDFluxJacobian(Q_barra(i+1,j_sol,1),Q_barra(i+1,j_sol,2),Q_barra(i+1,j_sol,3), &
                    Q_barra(i+1,j_sol,4),ksi_x(i+1,j_sol),ksi_y(i+1,j_sol),10.0d0**-5,A_barra)
            do index_i = 1, dim
                do index_j = 1, dim
                    A_plus(index_i,index_j,i-1) = 0.5d0*delta_t(i,j_sol)*A_barra(index_i,index_j) + L_ksi*Identy(index_i,index_j)
                end do
            end do
            
            L_ksi = dis_imp_ksi(i,j_sol,eps_dis_i,1)
            u = Q_barra(i-1,j_sol,2)/Q_barra(i-1,j_sol,1)
            v = Q_barra(i-1,j_sol,3)/Q_barra(i-1,j_sol,1)
call jacobian_ksi(u,v,Q_barra(i-1,j_sol,4),Q_barra(i-1,j_sol,1),ksi_x(i-1,j_sol),ksi_y(i-1,j_sol),dim,A_barra)
call FDFluxJacobian(Q_barra(i-1,j_sol,1),Q_barra(i-1,j_sol,2),Q_barra(i-1,j_sol,3), &
                    Q_barra(i-1,j_sol,4),ksi_x(i-1,j_sol),ksi_y(i-1,j_sol),10.0d0**-5,A_barra)            
            do index_i = 1, dim
                do index_j = 1, dim
                    A_minus(index_i,index_j,i-1) = -0.5d0*delta_t(i,j_sol)*A_barra(index_i,index_j) + L_ksi*Identy(index_i,index_j)
                end do
            end do
            
    end do

    ! now create the vector B_sys, i.e., vector b of the system
    !
    ! *********
    ! i really can dropp the j index because the subroutine call 
    ! is made for each value of the index j
    ! *********

        do i = 2, imax - 1
            call compute_residue(i,j_sol)
            Bx_sys(1,i-1) = -residue(i,j_sol,1)
            Bx_sys(2,i-1) = -residue(i,j_sol,2)
            Bx_sys(3,i-1) = -residue(i,j_sol,3)
            Bx_sys(4,i-1) = -residue(i,j_sol,4)
        end do
    
    ! we need to solve the block tridiagonal system j times
    ! blktriad(maind,lower,upper,id,md,xb,x)
    ! put artificial dissipation
    
    do i = 2, imax - 1
        L_ksi = dis_imp_ksi(i,j_sol,eps_dis_i,2)
        do index_i = 1, dim 
            do index_j = 1, dim
                main_x(index_i,index_j,i-1) = Id_x(index_i,index_j,i) + L_ksi*Identy(index_i,index_j)
            end do 
        end do 
    end do

    call blktriad(main_x,A_minus,A_plus,dim,imax-2,Bx_sys,deltaQ_til_i) 
    
    ! update j_sol
    
    do i = 2, imax - 1
        deltaQ_til(1,i,j_sol) = deltaQ_til_i(1,i-1)
        deltaQ_til(2,i,j_sol) = deltaQ_til_i(2,i-1)
        deltaQ_til(3,i,j_sol) = deltaQ_til_i(3,i-1)
        deltaQ_til(4,i,j_sol) = deltaQ_til_i(4,i-1)
    end do 
    ! j_sol = j_sol + 1

end do looping_j_sol
    
    ! step ii)
    
    !do while ( i_sol <= imax - 1 )
    do i_sol = 2, imax - 1
        
        do j = 2, jmax - 1
            
            u = Q_barra(i_sol,j+1,2)/Q_barra(i_sol,j+1,1)
            v = Q_barra(i_sol,j+1,3)/Q_barra(i_sol,j+1,1)
call jacobian_eta(u,v,Q_barra(i_sol,j+1,4),Q_barra(i_sol,j+1,1),eta_x(i_sol,j+1),eta_y(i_sol,j+1),dim,B_barra)
call FDFluxJacobian(Q_barra(i_sol,j+1,1),Q_barra(i_sol,j+1,2),Q_barra(i_sol,j+1,3), &
                    Q_barra(i_sol,j+1,4),eta_x(i_sol,j+1),eta_y(i_sol,j+1),10.0d0**-5,A_barra)            
            L_eta = dis_imp_eta(i_sol,j,eps_dis_i,3)

            do index_i = 1, dim
                do index_j = 1, dim
                    B_plus(index_i,index_j,j-1) = 0.5d0*delta_t(i_sol,j)*B_barra(index_i,index_j) + L_eta*Identy(index_i,index_j)
                end do
            end do
            
            u = Q_barra(i_sol,j-1,2)/Q_barra(i_sol,j-1,1)
            v = Q_barra(i_sol,j-1,3)/Q_barra(i_sol,j-1,1)
call jacobian_eta(u,v,Q_barra(i_sol,j-1,4),Q_barra(i_sol,j-1,1),eta_x(i_sol,j-1),eta_y(i_sol,j-1),dim,B_barra)
call FDFluxJacobian(Q_barra(i_sol,j-1,1),Q_barra(i_sol,j-1,2),Q_barra(i_sol,j-1,3), &
                    Q_barra(i_sol,j-1,4),eta_x(i_sol,j-1),eta_y(i_sol,j-1),10.0d0**-5,A_barra)            
            L_eta = dis_imp_eta(i_sol,j,eps_dis_i,1)
            
            do index_i = 1, dim
                do index_j = 1, dim
                    B_minus(index_i,index_j,j-1) = -0.5d0*delta_t(i_sol,j)*B_barra(index_i,index_j) + L_eta*Identy(index_i,index_j)
                end do
            end do
            
            
        end do
        
        ! now create By_sys
        
        do j = 2, jmax - 1
            By_sys(1,j-1) = deltaQ_til(1,i_sol,j)
            By_sys(2,j-1) = deltaQ_til(2,i_sol,j)
            By_sys(3,j-1) = deltaQ_til(3,i_sol,j)
            By_sys(4,j-1) = deltaQ_til(4,i_sol,j)
        end do
        
        ! solve the block tridiagonal i times

        do j = 2, jmax - 1
            L_eta = dis_imp_eta(i_sol,j,eps_dis_i,2)
            do index_i = 1, dim 
                do index_j = 1, dim
                    main_y(index_i,index_j,j-1) = Id_y(index_i,index_j,j) + L_eta*Identy(index_i,index_j)
                end do 
            end do 
        end do

        call blktriad(main_y,B_minus,B_plus,dim,jmax-2,By_sys,deltaQ_til_j)
        
        ! add one to the index loop
        
        do j = 2, jmax - 1
            deltaQ(i_sol,j,1) = deltaQ_til_j(1,j-1)        
            deltaQ(i_sol,j,2) = deltaQ_til_j(2,j-1)
            deltaQ(i_sol,j,3) = deltaQ_til_j(3,j-1)
            deltaQ(i_sol,j,4) = deltaQ_til_j(4,j-1)          
        end do 
        
        ! update Q_barra
        
        do j = 2, jmax - 1
            Q_barra(i_sol,j,1) = deltaQ(i_sol,j,1) + Q_barra(i_sol,j,1)
            Q_barra(i_sol,j,2) = deltaQ(i_sol,j,2) + Q_barra(i_sol,j,2)
            Q_barra(i_sol,j,3) = deltaQ(i_sol,j,3) + Q_barra(i_sol,j,3)
            Q_barra(i_sol,j,4) = deltaQ(i_sol,j,4) + Q_barra(i_sol,j,4)   
        end do 
        ! i_sol = i_sol + 1
    end do 

deallocate(Id_x)
deallocate(Bx_sys)
deallocate(By_sys)
!
!
!
deallocate(A_barra, B_barra, M_barra)
deallocate(A_minus, A_plus)
deallocate(B_minus, B_plus)
deallocate(deltaQ_til)
deallocate(deltaQ_til_i,deltaQ_til_j, deltaQ)
deallocate(Id_y)
deallocate(main_x,main_y)
deallocate(Identy)
!
!
!
!209 format (500f12.6)
end subroutine implicit_beam_warming
