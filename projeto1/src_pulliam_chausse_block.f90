subroutine pulliam_chausse_block
    use vars
    use diagonalization
    use functions
    implicit none
    real(8),dimension(:,:,:),allocatable      :: main,upper,lower
    real(8),dimension(:,:,:),allocatable      :: x1_f, x2_f, delta_Q
    real(8),dimension(:,:),allocatable        :: u,v,rho,p,a 
    real(8),dimension(:,:),allocatable        :: inv_t_xi
    real(8),dimension(:,:),allocatable        :: right_side
    real(8),dimension(:,:),allocatable        :: x1,x2
    real(8),dimension(:,:),allocatable        :: Identy
    real(8),dimension(:,:),allocatable        :: n_inverse, Teta
    real(8),dimension(:),allocatable          :: diag_plus, diag_minus
    real(8),dimension(:),allocatable          :: residue_pc, result
    real(8),dimension(:),allocatable          :: aux_mult
    real(8)                                   :: L_ksi, L_eta
    integer(4)                                :: index

allocate(u(imax,jmax),v(imax,jmax),rho(imax,jmax))
allocate(p(imax,jmax),a(imax,jmax))

allocate(aux_mult(dim),result(dim))
allocate(Identy(dim,dim))

do j = 1, jmax
        do i = 1, imax
            Q_dis(i,j,1) = Q_barra(i,j,1)/metric_jacobian(i,j)
            Q_dis(i,j,2) = Q_barra(i,j,2)/metric_jacobian(i,j)
            Q_dis(i,j,3) = Q_barra(i,j,3)/metric_jacobian(i,j)
            Q_dis(i,j,4) = Q_barra(i,j,4)/metric_jacobian(i,j)

            u(i,j)   = Q_barra(i,j,2)/Q_barra(i,j,1)
            v(i,j)   = Q_barra(i,j,3)/Q_barra(i,j,1)
            rho(i,j) = Q_barra(i,j,1)/metric_jacobian(i,j)
            p(i,j)   = (gama-1.0d0)*(Q_barra(i,j,4)/metric_jacobian(i,j) &
                       - 0.50d0*rho(i,j)*(u(i,j)**2.0d0+v(i,j)**2.0d0))
            a(i,j)   = sqrt(gama*p(i,j)/rho(i,j))
        end do
end do

! start to solve the system at the ksi direction

allocate(inv_t_xi(dim,dim))
allocate(diag_plus(dim),diag_minus(dim))
allocate(main(dim,dim,imax-2),upper(dim,dim,imax-2),lower(dim,dim,imax-2))
allocate(right_side(dim,imax-2))
allocate(x1(dim,imax-2))
allocate(x1_f(imax,jmax,dim))
allocate(residue_pc(dim))

main  = 0.0d0
upper = 0.0d0
lower = 0.0d0
right_side = 0.0d0
x1 = 0.0d0
inv_t_xi = 0.0d0

do j = 2, jmax - 1

    do i = 2, imax - 1 

        inv_t_xi           = inv_T_ksi(u(i,j),v(i,j),rho(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
        
        call compute_residue(i,j)
        residue_pc(1:dim)     = -residue(i,j,1:dim)
        aux_mult              = matmul(inv_T_xi,residue_pc)
        right_side(1:dim,i-1) = aux_mult(1:dim)

        diag_plus  = diag_ksi(U_contravariant(i+1,j),a(i+1,j),ksi_x(i+1,j),ksi_y(i+1,j),dim)
        
        L_ksi = dis_imp_ksi(i,j,eps_dis_i,3)
        do index = 1, dim
            upper(index,index,i-1)  =  0.50d0*delta_t(i,j)*diag_plus(index) + L_ksi*Identy(index,index)
        end do

        diag_minus = diag_ksi(U_contravariant(i-1,j),a(i-1,j),ksi_x(i-1,j),ksi_y(i-1,j),dim)

        L_ksi = dis_imp_ksi(i,j,eps_dis_i,1)
        do index = 1, dim
            lower(index,index,i-1) = -0.50d0*delta_t(i,j)*diag_minus(index) + L_ksi*Identy(index,index)
        end do

        L_ksi = dis_imp_ksi(i,j,eps_dis_i,2)
        do index = 1, dim
            main(index,index,i-1) = 1.0d0 + L_ksi
        end do

    end do 

    call blktriad(main,lower,upper,dim,imax-2,right_side,x1) 

    do i = 2, imax - 1 
        x1_f(i,j,1:dim) = x1(1:dim,i-1)
    end do 

end do 

deallocate(inv_t_xi)
deallocate(diag_plus,diag_minus)
deallocate(main,upper,lower)
deallocate(right_side)
deallocate(residue_pc)

! start to solve the system at the eta-direction

allocate(diag_plus(dim),diag_minus(dim))
allocate(main(dim,dim,jmax-2),upper(dim,dim,jmax-2),lower(dim,dim,jmax-2))
allocate(right_side(dim,jmax-2))
allocate(x2(dim,jmax-2))
allocate(x2_f(imax,jmax,dim))
allocate(n_inverse(dim,dim))

main  = 0.0d0
upper = 0.0d0
lower = 0.0d0
right_side = 0.0d0
x2 = 0.0d0
n_inverse = 0.0d0

do i = 2, imax - 1

    do j = 2, jmax - 1

        aux_mult(1:dim) = x1_f(i,j,1:dim)
        n_inverse = inv_N_matrix(ksi_x(i,j),ksi_y(i,j),eta_x(i,j),eta_y(i,j),dim)
        result = matmul(n_inverse,aux_mult)
        right_side(1:dim,j-1) = result(1:dim)

        diag_plus = diag_eta(V_contravariant(i,j+1),a(i,j+1),eta_x(i,j+1),eta_y(i,j+1),dim)

        L_eta = dis_imp_eta(i,j,eps_dis_i,3)
        do index = 1, dim
            upper(index,index,j-1)  =  0.50d0*delta_t(i,j)*diag_plus(index) + L_eta*Identy(index,index)
        end do

        diag_minus = diag_eta(V_contravariant(i,j-1),a(i,j-1),eta_x(i,j-1),eta_y(i,j-1),dim)

        L_eta = dis_imp_eta(i,j,eps_dis_i,1)
        do index = 1, dim
            lower(index,index,j-1) = -0.50d0*delta_t(i,j)*diag_minus(index) + L_eta*Identy(index,index)
        end do

        L_eta = dis_imp_eta(i,j,eps_dis_i,2)
        do index = 1, dim
            main(index,index,j-1) = 1.0d0 + L_eta
        end do

    end do 

    call blktriad(main,lower,upper,dim,jmax-2,right_side,x2) 

    do j = 2, jmax - 1 
        x2_f(i,j,1:dim) = x2(1:dim,j-1)
    end do 

end do

deallocate(diag_plus,diag_minus)
deallocate(main,upper,lower)
deallocate(right_side)
deallocate(x1)
deallocate(x2)

deallocate(Identy)

! calculate delta_Q in order to update the numerical solution

allocate(Teta(dim,dim))
allocate(delta_Q(imax,jmax,dim))

Teta = 0.0d0
delta_Q = 0.0d0

do j = 2, jmax - 1
    do i = 2, imax - 1
        Teta = T_eta(u(i,j),v(i,j),rho(i,j),a(i,j),eta_x(i,j),eta_y(i,j),dim)
        aux_mult(1:dim) = x2_f(i,j,1:dim)
        result = matmul(Teta,aux_mult)
        delta_Q(i,j,1:dim) = result(1:dim)
    end do 
end do 

deallocate(u,v,rho)
deallocate(p,a)
deallocate(x2_f)
deallocate(x1_f)
deallocate(n_inverse,Teta)

! update the solution

do j = 2, jmax - 1
    do i = 2, imax - 1

        Q_barra(i,j,1) = Q_barra(i,j,1) + delta_Q(i,j,1) 
        Q_barra(i,j,2) = Q_barra(i,j,2) + delta_Q(i,j,2) 
        Q_barra(i,j,3) = Q_barra(i,j,3) + delta_Q(i,j,3) 
        Q_barra(i,j,4) = Q_barra(i,j,4) + delta_Q(i,j,4) 
        
    end do
end do

deallocate(delta_Q)
deallocate(aux_mult,result)

end subroutine pulliam_chausse_block