subroutine pulliam_chausse
    use vars
    use diagonalization
    use functions
    implicit none
    !integer(4)                                     :: i_sol, j_sol
    real(8)                                        :: rho_t
    real(8),dimension(:,:),allocatable             :: inverse_T_ksi, T_eta_matrix
    real(8),dimension(:,:,:),allocatable           :: x1_pc, x2_pc
    real(8),dimension(:,:),allocatable             :: N_pc
    real(8),dimension(:),allocatable               :: aux_residue, aux_x2
    real(8),dimension(:),allocatable               :: result, ans
    real(8),dimension(:),allocatable               :: diag_minus,diag_plus
    real(8),dimension(:),allocatable               :: deltaQ
    reaL(8),dimension(:),allocatable               :: lower, main, upper, d_thomas
    real(8)                                        :: L_ksi, L_eta
    integer(4)                                     :: ii, jj, index_i, index_j, index
!
!
allocate(inverse_T_ksi(dim,dim),T_eta_matrix(dim,dim))
allocate(aux_residue(dim),x1_pc(imax,jmax,dim))
allocate(x2_pc(imax,jmax,dim))
allocate(result(dim),N_pc(dim,dim))
allocate(diag_minus(dim),diag_plus(dim))
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
allocate(lower(dim*(imax-2)),main(dim*(imax-2)),upper(dim*(imax-2)))
allocate(d_thomas(dim*(imax-2)),ans(dim*(imax-2)))
!
! system for ksi 
!
do j = 2, jmax - 1
    do i = 2, imax - 1
        !
        !
            rho_t         = Q_barra(i,j,1)/metric_jacobian(i,j)
            inverse_T_ksi = inv_T_ksi(u(i,j),v(i,j),rho_t,a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
                call compute_residue(i,j)
                aux_residue(1) = -residue(i,j,1)
                aux_residue(2) = -residue(i,j,2)
                aux_residue(3) = -residue(i,j,3)
                aux_residue(4) = -residue(i,j,4)
            result = matmul(inverse_T_ksi,aux_residue)
            diag_plus  = diag_ksi(U_contravariant(i+1,j),a(i+1,j),ksi_x(i+1,j),ksi_y(i+1,j),dim)
            diag_minus = diag_ksi(U_contravariant(i-1,j),a(i-1,j),ksi_x(i-1,j),ksi_y(i-1,j),dim)
            !    a - sub-diagonal (means it is the diagonal below the main diagonal)
            !    b - the main diagonal
            !    c - sup-diagonal (means it is the diagonal above the main diagonal)
            !    d - right part
            !    x - the answer
            !    n - number of equations
            L_ksi = dis_imp_ksi(i,j,eps_dis_i,3)
            diag_plus  = 0.50d0*delta_t(i,j)*diag_plus + L_ksi
            !
            L_ksi = dis_imp_ksi(i,j,eps_dis_i,1)
            diag_minus = -0.50d0*delta_t(i,j)*diag_minus + L_ksi
            !
            L_ksi = dis_imp_ksi(i,j,eps_dis_i,2)
            do index = 1, dim 
                    ii = index + (i-2)*dim
                    lower(ii)    = diag_minus(index)
                    main(ii)     = 1.0d0 + L_ksi
                    upper(ii)    = diag_plus(index)
                    d_thomas(ii) = result(index)  
            end do
            !
            !
    end do 
        call thomas_pulliam_chausse(lower,main,upper,d_thomas,ans,dim*(imax-2))
        do index_i = 2, imax - 1
            do index = 1, dim
                x1_pc(index_i,j,index) = ans(index +(index_i-2)*dim)
            end do
        end do
end do
!
!
deallocate(lower,main,upper)
deallocate(d_thomas,ans)
!
!
allocate(lower(dim*(jmax-2)),main(dim*(jmax-2)),upper(dim*(jmax-2)))
allocate(d_thomas(dim*(jmax-2)),ans(dim*(jmax-2)))
!
!
do i = 2, imax - 1
    do j = 2, jmax - 1
        !
        !
        N_pc          = inv_N_matrix(ksi_x(i,j),ksi_y(i,j),eta_x(i,j),eta_y(i,j),dim)
        !
        !
        aux_residue(1) = x1_pc(i,j,1) 
        aux_residue(2) = x1_pc(i,j,2)
        aux_residue(3) = x1_pc(i,j,3)
        aux_residue(4) = x1_pc(i,j,4)
        result = matmul(N_pc,aux_residue)
        !
        !
        diag_plus  = diag_eta(V_contravariant(i,j+1),a(i,j+1),eta_x(i,j+1),eta_y(i,j+1),dim)
        diag_minus = diag_eta(V_contravariant(i,j-1),a(i,j-1),eta_x(i,j-1),eta_y(i,j-1),dim)
        !
        !
        L_eta = dis_imp_eta(i,j,eps_dis_i,3)
        diag_plus  = 0.50d0*delta_t(i,j)*diag_plus + L_eta 
        !
        L_eta = dis_imp_eta(i,j,eps_dis_i,1)
        diag_minus = -0.50d0*delta_t(i,j)*diag_minus + L_eta
        !
        L_eta = dis_imp_eta(i,j,eps_dis_i,2)
        do index = 1, dim 
                jj = index + (j-2)*dim 
                lower(jj)    = diag_minus(index)
                main(jj)     = 1.0d0 + L_eta
                upper(jj)    = diag_plus(index)
                d_thomas(jj) = result(index)  
        end do
        !
        !
    end do 
        call thomas_pulliam_chausse(lower,main,upper,d_thomas,ans,dim*(jmax-2))   
        do index_j = 2, jmax - 1
            do index = 1, dim
                x2_pc(i,index_j,index) = ans(index + (index_j-2)*dim)
            end do
        end do
end do 
!
!
deallocate(lower,main,upper)
deallocate(d_thomas)
deallocate(inverse_T_ksi)
deallocate(aux_residue,x1_pc)
deallocate(result,N_pc)
deallocate(diag_minus,diag_plus)
!
! update the q
!
allocate(deltaQ(dim),aux_x2(dim))
!
!   
do j = 2, jmax - 1
    do i = 2, imax - 1
        !
        aux_x2(1) = x2_pc(i,j,1)
        aux_x2(2) = x2_pc(i,j,2)
        aux_x2(3) = x2_pc(i,j,3)
        aux_x2(4) = x2_pc(i,j,4)
        !
        rho_t         = Q_barra(i,j,1)/metric_jacobian(i,j)
        T_eta_matrix  = T_eta(u(i,j),v(i,j),rho_t,a(i,j),eta_x(i,j),eta_y(i,j),dim)
        deltaQ = matmul(T_eta_matrix,aux_x2)
        !
        Q_barra(i,j,1) = deltaQ(1) + Q_barra(i,j,1)         
        Q_barra(i,j,2) = deltaQ(2) + Q_barra(i,j,2)
        Q_barra(i,j,3) = deltaQ(3) + Q_barra(i,j,3)
        Q_barra(i,j,4) = deltaQ(4) + Q_barra(i,j,4)
        !
    end do 
end do 
!
!
deallocate(x2_pc,deltaQ,aux_x2,T_eta_matrix)
!
!
end subroutine pulliam_chausse