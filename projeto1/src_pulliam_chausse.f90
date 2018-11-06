subroutine pulliam_chausse
    use vars
    use diagonalization
    use functions
    implicit none
    real(8),dimension(:),allocatable                   :: diag_plus, diag_minus
    real(8),dimension(:),allocatable                   :: result, aux_mult
    real(8),dimension(:),allocatable                   :: lower,main,upper
    real(8),dimension(:),allocatable                   :: x_sys,d_sys
    real(8),dimension(:,:),allocatable                 :: inv_t_xi , n_inverse
    real(8),dimension(:,:),allocatable                 :: Teta
    real(8),dimension(:,:,:),allocatable               :: x1_sys, right_side
    real(8)                                            :: rho_t, L_ksi, L_eta
    integer(4)                                         :: index, i_sol, j_sol
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
!
allocate(x1_sys(imax,jmax,dim))
allocate(right_side(imax,jmax,dim))
allocate(inv_t_xi(dim,dim),n_inverse(dim,dim))
allocate(Teta(dim,dim))
allocate(result(dim))
allocate(diag_minus(dim),diag_plus(dim))
!
!
allocate(lower(imax-2),main(imax-2),upper(imax-2))
allocate(x_sys(imax-2),d_sys(imax-2),aux_mult(dim))
!
!
diag_minus = 0.0d0
diag_plus  = 0.0d0
d_sys = 0.0d0
x1_sys = 0.0d0
x_sys  = 0.0d0
lower  = 0.0d0
upper  = 0.0d0
main   = 1.0d0
aux_mult = 0.0d0
right_side = 0.0d0
!
!
do j = 2, jmax - 1
    do i = 2, imax - 1
        call compute_residue(i,j)
        aux_mult(1:dim) = -residue(i,j,1:dim)
        rho_t = Q_barra(i,j,1)/metric_jacobian(i,j)
        inv_t_xi = inv_T_ksi(u(i,j),v(i,j),rho_t,a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
        result = matmul(inv_t_xi,aux_mult)
        right_side(i,j,1:dim) = result(1:dim)
    end do 
end do
!
do index = 1, dim
    !
    do j_sol = 2, jmax - 1
        !
        do i = 2, imax - 1
            !
            d_sys(i-1) = right_side(i,j_sol,index)
            !
            L_ksi = dis_imp_ksi(i,j_sol,eps_dis_i,3)
            diag_plus  = diag_ksi(U_contravariant(i+1,j_sol),a(i+1,j_sol),ksi_x(i+1,j_sol),ksi_y(i+1,j_sol),dim)
                upper(i-1) = 0.50d0*delta_t(i,j_sol)*diag_plus(index) + L_ksi
            L_ksi = dis_imp_ksi(i,j_sol,eps_dis_i,1)
            diag_minus = diag_ksi(U_contravariant(i-1,j_sol),a(i-1,j_sol),ksi_x(i-1,j_sol),ksi_y(i-1,j_sol),dim)
                lower(i-1) = -0.50d0*delta_t(i,j_sol)*diag_minus(index) + L_ksi
            L_ksi = dis_imp_ksi(i,j_sol,eps_dis_i,2)
                main(i-1) = main(i-1) + L_ksi
            !
        end do 
        !    a - sub-diagonal (means it is the diagonal below the main diagonal)
        !    b - the main diagonal
        !    c - sup-diagonal (means it is the diagonal above the main diagonal)
        !    d - right part
        !    x - the answer
        !    n - number of equations
        call thomas_pulliam_chausse(lower,main,upper,d_sys,x_sys,imax-2)
        !
        do i = 2, imax - 1
            x1_sys(i,j_sol,index) = x_sys(i-1)
        end do
        !
    end do
    !
end do
!
!
deallocate(lower,main,upper)
deallocate(x_sys,d_sys)
!
! start to solve the system in the eta location
!
allocate(lower(jmax-2),main(jmax-2),upper(jmax-2))
allocate(x_sys(jmax-2),d_sys(jmax-2))
!
!
d_sys = 0.0d0
x_sys  = 0.0d0
lower  = 0.0d0
upper  = 0.0d0
main   = 1.0d0
aux_mult = 0.0d0
right_side = 0.0d0
!
!
do j = 2, jmax - 1
    do i = 2, imax - 1
        aux_mult(1:dim) = x1_sys(i,j,1:dim)
        n_inverse = inv_N_matrix(ksi_x(i,j),ksi_y(i,j),eta_x(i,j),eta_y(i,j),dim)
        result = matmul(n_inverse,aux_mult)
        right_side(i,j,1:dim) = result(1:dim)
    end do 
end do
!
do index = 1, dim
    ! do i = 2, imax - 1
    do i_sol = 2, imax - 1
        !
        do j = 2, jmax - 1
            !
            d_sys(j-1) = right_side(i_sol,j,index)
            !
            L_eta = dis_imp_eta(i_sol,j,eps_dis_i,3)
            diag_plus = diag_eta(V_contravariant(i_sol,j+1),a(i_sol,j+1),eta_x(i_sol,j+1),eta_y(i_sol,j+1),dim)
                upper(j-1) = 0.50d0*delta_t(i_sol,j)*diag_plus(index) + L_eta
            diag_minus = diag_eta(V_contravariant(i_sol,j-1),a(i_sol,j-1),eta_x(i_sol,j-1),eta_y(i_sol,j-1),dim)
            L_eta = dis_imp_eta(i_sol,j,eps_dis_i,1)
                lower(j-1) = -0.50d0*delta_t(i_sol,j)*diag_minus(index) + L_eta
            L_eta = dis_imp_eta(i_sol,j,eps_dis_i,2)
                main(j-1) = main(j-1) + L_eta
            !
        end do
        !
        call thomas_pulliam_chausse(lower,main,upper,d_sys,x_sys,jmax-2)
        do j = 2, jmax - 1
            x1_sys(i_sol,j,index) = x_sys(j-1)
            ! if (index == 1) write(*,*) x1_sys(i,2,index), x_sys(1)
        end do
        !
    end do
    !
end do 
!
!
deallocate(lower,main,upper)
deallocate(x_sys,d_sys)
deallocate(diag_minus,diag_plus)
!
!
do j = 2, jmax - 1
    do i = 2, imax - 1
        rho_t = Q_barra(i,j,1)/metric_jacobian(i,j)
        Teta = T_eta(u(i,j),v(i,j),rho_t,a(i,j),eta_x(i,j),eta_y(i,j),dim)
        aux_mult(1:dim) = x1_sys(i,j,1:dim)
        result = matmul(Teta,aux_mult)
        right_side(i,j,1:dim) = result(1:dim)
    end do 
end do 
!
deallocate(x1_sys)
deallocate(inv_t_xi,n_inverse)
deallocate(result,aux_mult,Teta)
!
! update the solution vector
!
do j = 2, jmax - 1
    do i = 2, imax - 1
        !
        Q_barra(i,j,1) = Q_barra(i,j,1) + right_side(i,j,1) 
        Q_barra(i,j,2) = Q_barra(i,j,2) + right_side(i,j,2) 
        Q_barra(i,j,3) = Q_barra(i,j,3) + right_side(i,j,3) 
        Q_barra(i,j,4) = Q_barra(i,j,4) + right_side(i,j,4) 
        !
    end do
end do
!
!
        if (iter == 50) then 
            open(999,file='pc')
            do j = 2, jmax - 1
                do i = 2, imax - 1 
                    write(999,*) iter,i,j,Q_barra(i,j,1)/metric_jacobian(i,j)     
                end do 
            end do
            close(999)    
        end if
!
deallocate(right_side)
!
!
end subroutine pulliam_chausse