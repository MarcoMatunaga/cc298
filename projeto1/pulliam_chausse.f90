subroutine pulliam_chausse
    use vars
    use diagonalization
    implicit none
    !integer(4)                                     :: i_sol, j_sol
    real(8)                                        :: rho_t
    real(8),dimension(:,:),allocatable             :: inverse_T_ksi, T_eta_matrix
    real(8),dimension(:,:,:),allocatable           :: x1_pc, x2_pc
    real(8),dimension(:,:),allocatable             :: lambda_plus, lambda_minus, N_pc
    real(8),dimension(:),allocatable               :: Identy
    real(8),dimension(:),allocatable               :: aux_residue, aux_x2
    real(8),dimension(:),allocatable               :: result, ans
    real(8),dimension(:),allocatable               :: diag_minus,diag_plus
    real(8),dimension(:),allocatable               :: deltaQ
!
!
allocate(Identy(dim))
allocate(inverse_T_ksi(dim,dim),T_eta_matrix(dim,dim))
allocate(lambda_minus(dim,dim),lambda_plus(dim,dim))
allocate(aux_residue(dim),x1_pc(imax,jmax,dim))
allocate(x2_pc(imax,jmax,dim))
allocate(result(dim),ans(dim),N_pc(dim,dim))
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
do i = 1, dim
     Identy(i) = 1.0d0 
end do
!
!
do j = 2, jmax - 1
    !
    !
    do i = 2, imax - 1
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
        diag_plus  = 0.50d0*delta_t(i,j)*diag_plus
        diag_minus = -0.50d0*delta_t(i,j)*diag_minus
        call thomas_pulliam_chausse(diag_minus,Identy,diag_plus,result,ans,dim)
        !
        !
        x1_pc(i,j,1) = ans(1)
        x1_pc(i,j,2) = ans(2)
        x1_pc(i,j,3) = ans(3)
        x1_pc(i,j,4) = ans(4)
        !
        !
    end do
    !
    !
end do
!
!
do i = 2, imax -1
    !
    !
    do j = 2, jmax - 1
        N_pc          = inv_N_matrix(ksi_x(i,j),ksi_y(i,j),eta_x(i,j),eta_y(i,j),dim)
        !
        !
        aux_residue(1) = x1_pc(i,j,1)
        aux_residue(2) = x1_pc(i,j,2)
        aux_residue(3) = x1_pc(i,j,3)
        aux_residue(4) = x1_pc(i,j,4)
        !
        !
        result = matmul(N_pc,aux_residue)
        !
        !
        diag_plus  = diag_eta(V_contravariant(i,j+1),a(i,j+1),eta_x(i,j+1),eta_y(i,j+1),dim)
        diag_minus = diag_eta(V_contravariant(i,j-1),a(i,j-1),eta_x(i,j-1),eta_y(i,j-1),dim)
        !
        !
        diag_plus  = 0.50d0*delta_t(i,j)*diag_plus
        diag_minus = -0.50d0*delta_t(i,j)*diag_minus
        call thomas_pulliam_chausse(diag_minus,Identy,diag_plus,result,ans,dim)
        !
        !
        x2_pc(i,j,1) = ans(1)
        x2_pc(i,j,2) = ans(2)
        x2_pc(i,j,3) = ans(3)
        x2_pc(i,j,4) = ans(4)
        !
        !
    end do
    !
    !
end do
!
!
deallocate(Identy)
deallocate(inverse_T_ksi)
deallocate(lambda_minus,lambda_plus)
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