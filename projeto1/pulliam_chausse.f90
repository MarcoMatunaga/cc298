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
    real(8)                                            :: rho_t
    integer(4)                                         :: index
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
allocate(lower(imax),main(imax),upper(imax))
allocate(x_sys(imax),d_sys(imax),aux_mult(dim))
!
!
diag_minus = 0.0d0
diag_plus  = 0.0d0
d_sys = 0.0d0
x1_sys = 0.0d0
x_sys  = 0.0d0
lower  = 0.0d0
upper  = 0.0d0
main   = 0.0d0
aux_mult = 0.0d0
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
index = 1
!
do while (index <= dim)
        !
        do j = 2, jmax - 1
        !
                i = 1
                !
                diag_plus  = diag_ksi(U_contravariant(i+1,j),a(i+1,j),ksi_x(i+1,j),ksi_y(i+1,j),dim)
                       upper(i) = 0.50d0*delta_t(i,j)*diag_plus(index)
                       main(i)  = 1.0d0
                !
            do i = 2, imax - 1
                !
                        d_sys(i) = right_side(i,j,index)
                !
                diag_plus  = diag_ksi(U_contravariant(i+1,j),a(i+1,j),ksi_x(i+1,j),ksi_y(i+1,j),dim)
                        upper(i) = 0.50d0*delta_t(i,j)*diag_plus(index)
                diag_minus = diag_ksi(U_contravariant(i-1,j),a(i-1,j),ksi_x(i-1,j),ksi_y(i-1,j),dim)
                        lower(i) = -0.50d0*delta_t(i,j)*diag_minus(index)
                !
                        main(i) = 1.0d0
            end do 
                !
                i = imax
                !
                diag_minus = diag_ksi(U_contravariant(i-1,j),a(i-1,j),ksi_x(i-1,j),ksi_y(i-1,j),dim)
                        lower(i) = -0.50d0*delta_t(i,j)*diag_plus(index)
                        main(i)  = 1.0d0
                !
            !    a - sub-diagonal (means it is the diagonal below the main diagonal)
            !    b - the main diagonal
            !    c - sup-diagonal (means it is the diagonal above the main diagonal)
            !    d - right part
            !    x - the answer
            !    n - number of equations
            call thomas_pulliam_chausse(lower,main,upper,d_sys,x_sys,imax)
            do i = 2, imax - 1
                    !write(*,*) iter,i,j,lower(index+(i-1)*dim),main(index+(i-1)*dim),&
                    !          upper(index+(i-1)*dim),x_sys(index+(i-1)*dim)
                    x1_sys(i,j,index) = x_sys(i)
            end do
        !
        end do
        !
    index = index + 1
end do
!
!
deallocate(lower,main,upper)
deallocate(x_sys,d_sys)
!
! start to solve the system in the eta location
!
allocate(lower(jmax),main(jmax),upper(jmax))
allocate(x_sys(jmax),d_sys(jmax))
!
!
d_sys = 0.0d0
x_sys  = 0.0d0
lower  = 0.0d0
upper  = 0.0d0
main   = 0.0d0
aux_mult = 0.0d0
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
index = 1
!
do while (index <= dim)
    do i = 2, imax - 1
        !
        j = 1
        !
        diag_plus  = diag_eta(V_contravariant(i,j+1),a(i,j+1),eta_x(i,j+1),eta_y(i,j+1),dim)
                upper(j) = 0.50d0*delta_t(i,j)*diag_plus(index)
                main(j)  = 1.0d0
        !
        do j = 2, jmax - 1
            !
            d_sys(j) = right_side(i,j,index)
            !
            diag_plus  = diag_eta(V_contravariant(i,j+1),a(i,j+1),eta_x(i,j+1),eta_y(i,j+1),dim)
                upper(j) = 0.50d0*delta_t(i,j)*diag_plus(index)
            diag_minus = diag_eta(V_contravariant(i,j-1),a(i,j-1),eta_x(i,j-1),eta_y(i,j-1),dim)
                lower(j) = -0.50d0*delta_t(i,j)*diag_minus(index)
                main(j) = 1.0d0
            !
        end do
        !
        j = jmax
        !
        diag_minus = diag_eta(V_contravariant(i,j-1),a(i,j-1),eta_x(i,j-1),eta_y(i,j-1),dim)
                lower(j) = -0.50d0*delta_t(i,j)*diag_plus(index)
                main(j)  = 1.0d0
        !
        call thomas_pulliam_chausse(lower,main,upper,d_sys,x_sys,jmax)
        do j = 2, jmax - 1
                !write(*,*) iter,i,j,lower(index+(j-1)*dim),main(index+(j-1)*dim),&
                !           upper(index+(j-1)*dim),x_sys(index+(j-1)*dim)
                x1_sys(i,j,index) = x_sys(j)
        end do
        !
    end do
    write(*,*) index
    index = index + 1
end do 
!
!
deallocate(lower,main,upper)
deallocate(x_sys,d_sys)
!
!
do j = 2, jmax - 1
    do i = 2, imax - 1
        rho_t = Q_barra(i,j,1)/metric_jacobian(i,j)
        Teta = T_eta(u(i,j),v(i,j),rho_t,a(i,j),eta_x(i,j),eta_y(i,j),dim)
        !
        do index = 1, dim
            aux_mult(index) = x1_sys(i,j,index)
        end do
        !
        result = matmul(Teta,aux_mult)
        !
        do index = 1, dim 
            Q_barra(i,j,index) = result(index) + Q_barra(i,j,index)
        end do
        !
    end do
end do
deallocate(x1_sys)
deallocate(inv_t_xi,n_inverse)
deallocate(result,aux_mult)
deallocate(diag_minus,diag_plus)
!
!
end subroutine pulliam_chausse