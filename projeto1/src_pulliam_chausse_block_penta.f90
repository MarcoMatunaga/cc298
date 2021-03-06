subroutine pulliam_chausse_block_penta
    use vars
    use diagonalization
    use functions
    implicit none
    real(8),dimension(:,:),allocatable        :: Diss_ksi_pc, Diss_eta_pc
    real(8),dimension(5)                      :: aux_diss
    real(8),dimension(:,:,:),allocatable      :: main,upper,lower
    real(8),dimension(:,:,:),allocatable      :: upper1,lower1
    real(8),dimension(:,:,:),allocatable      :: x1_f, x2_f
    real(8),dimension(:,:),allocatable        :: u,v,rho,p,a 
    real(8),dimension(:,:),allocatable        :: inv_t_xi
    real(8),dimension(:,:),allocatable        :: right_side
    real(8),dimension(:,:),allocatable        :: x1,x2
    real(8),dimension(:,:),allocatable        :: Identy
    real(8),dimension(:,:),allocatable        :: n_inverse, Teta
    real(8),dimension(:),allocatable          :: diag_plus, diag_minus
    real(8),dimension(:),allocatable          :: residue_pc, result
    real(8),dimension(:),allocatable          :: aux_mult
    integer(4)                                :: index

allocate(u(imax,jmax),v(imax,jmax),rho(imax,jmax))
allocate(p(imax,jmax),a(imax,jmax))
allocate(Diss_ksi_pc(5,dim),Diss_eta_pc(5,dim))

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
allocate(upper1(dim,dim,imax-2),lower1(dim,dim,imax-2))
allocate(right_side(dim,imax-2))
allocate(x1(dim,imax-2))
allocate(x1_f(imax,jmax,dim))
allocate(residue_pc(dim))

main   = 0.0d0
upper  = 0.0d0
lower  = 0.0d0
upper1 = 0.0d0
lower1 = 0.0d0
right_side = 0.0d0
x1 = 0.0d0
inv_t_xi = 0.0d0

do j = 2, jmax - 1

    do i = 2, imax - 1 

        !******* use functions to calculate the primitive variables and other variables
        inv_t_xi           = inv_T_ksi(u(i,j),v(i,j),rho(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
        
        call compute_residue(i,j)
        residue_pc(1:dim)     = -residue(i,j,1:dim)
        aux_mult              = matmul(inv_T_xi,residue_pc)
        right_side(1:dim,i-1) = aux_mult(1:dim)

        ! calculate the dissipation terms
            do index = 1, dim
                aux_diss = dis_ksi4_imp(i,j,index,eps_dis_e) 
                Diss_ksi_pc(1:5,index) = aux_diss(1:5)
            end do 

        do index = 1, dim
            upper1(index,index,i-1) = Diss_ksi_pc(1,index)
        end do 

        do index = 1, dim 
            lower1(index,index,i-1) = Diss_ksi_pc(5,index)
        end do 

        diag_plus  = diag_ksi(U_contravariant(i+1,j),a(i+1,j),ksi_x(i+1,j),ksi_y(i+1,j),dim)
        do index = 1, dim
            upper(index,index,i-1)  =  0.50d0*delta_t(i,j)*diag_plus(index) + Diss_ksi_pc(2,index)
        end do

        diag_minus = diag_ksi(U_contravariant(i-1,j),a(i-1,j),ksi_x(i-1,j),ksi_y(i-1,j),dim)
        do index = 1, dim
            lower(index,index,i-1) = -0.50d0*delta_t(i,j)*diag_minus(index) + Diss_ksi_pc(4,index)
        end do

        do index = 1, dim
            main(index,index,i-1) = 1.0d0 + Diss_ksi_pc(3,index)
        end do

    end do 

    call penta_block(lower1,lower,main,upper,upper1,right_side,dim,imax-2,x1)

    do i = 2, imax - 1 
        x1_f(i,j,1:dim) = x1(1:dim,i-1)
    end do 

end do 

deallocate(diag_plus,diag_minus)
deallocate(main,upper,lower)
deallocate(upper1,lower1)
deallocate(right_side)
deallocate(residue_pc)
deallocate(Diss_ksi_pc)

! start to solve the system at the eta-direction

allocate(diag_plus(dim),diag_minus(dim))
allocate(main(dim,dim,jmax-2),upper(dim,dim,jmax-2),lower(dim,dim,jmax-2))
allocate(upper1(dim,dim,jmax-2),lower1(dim,dim,jmax-2))
allocate(right_side(dim,jmax-2))
allocate(x2(dim,jmax-2))
allocate(x2_f(imax,jmax,dim))
allocate(n_inverse(dim,dim))

main   = 0.0d0
upper  = 0.0d0
lower  = 0.0d0
upper1 = 0.0d0
lower1 = 0.0d0
right_side = 0.0d0
x2 = 0.0d0
n_inverse = 0.0d0

do i = 2, imax - 1

    do j = 2, jmax - 1

        aux_mult(1:dim) = x1_f(i,j,1:dim)
        n_inverse = inv_N_matrix(ksi_x(i,j),ksi_y(i,j),eta_x(i,j),eta_y(i,j),dim)
        result = matmul(n_inverse,aux_mult) 

        right_side(1:dim,j-1) = result(1:dim)
        ! calculate the dissipation terms
            do index = 1, dim
                aux_diss = dis_eta4_imp(i,j,index,eps_dis_e) 
                Diss_eta_pc(1:5,index) = aux_diss(1:5)
            end do 

        do index = 1, dim
            upper1(index,index,j-1) = Diss_eta_pc(1,index)
        end do 

        do index = 1, dim 
            lower1(index,index,j-1) = Diss_eta_pc(5,index)
        end do 

        diag_plus = diag_eta(V_contravariant(i,j+1),a(i,j+1),eta_x(i,j+1),eta_y(i,j+1),dim)
        do index = 1, dim
            upper(index,index,j-1)  =  0.50d0*delta_t(i,j)*diag_plus(index) + Diss_eta_pc(2,index)
        end do

        diag_minus = diag_eta(V_contravariant(i,j-1),a(i,j-1),eta_x(i,j-1),eta_y(i,j-1),dim)
        do index = 1, dim
            lower(index,index,j-1) = -0.50d0*delta_t(i,j)*diag_minus(index) + Diss_eta_pc(4,index)
        end do

        do index = 1, dim
            main(index,index,j-1) = 1.0d0 + Diss_eta_pc(3,index)
        end do

    end do 

    call penta_block(lower1,lower,main,upper,upper1,right_side,dim,jmax-2,x2)

    do j = 2, jmax - 1 
        x2_f(i,j,1:dim) = x2(1:dim,j-1)
    end do 

end do

deallocate(diag_plus,diag_minus)
deallocate(main,upper,lower)
deallocate(upper1,lower1)
deallocate(right_side)
deallocate(x1)
deallocate(x2)
deallocate(Diss_eta_pc)

deallocate(Identy)

! calculate delta_Q in order to update the numerical solution

allocate(Teta(dim,dim))

Teta = 0.0d0
result = 0.0d0

do j = 2, jmax - 1
    do i = 2, imax - 1

        Teta = T_eta(u(i,j),v(i,j),rho(i,j),a(i,j),eta_x(i,j),eta_y(i,j),dim)
        aux_mult(1:dim) = x2_f(i,j,1:dim)
        result = matmul(Teta,aux_mult)

        ! update the solution

        Q_barra(i,j,1) = Q_barra(i,j,1) + result(1) 
        Q_barra(i,j,2) = Q_barra(i,j,2) + result(2) 
        Q_barra(i,j,3) = Q_barra(i,j,3) + result(3) 
        Q_barra(i,j,4) = Q_barra(i,j,4) + result(4) 

    end do 
end do 

deallocate(u,v,rho)
deallocate(p,a)
deallocate(x2_f)
deallocate(x1_f)
deallocate(n_inverse,Teta)
deallocate(aux_mult,result)
deallocate(inv_t_xi)

end subroutine pulliam_chausse_block_penta