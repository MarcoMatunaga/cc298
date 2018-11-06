subroutine sw_2nd
    use vars
    use diagonalization
    use fluxes_pos_neg
    implicit none
    real(8),dimension(:),allocatable             :: eig
    real(8),dimension(:,:),allocatable           :: diag_pos, diag_neg
    real(8),dimension(:,:,:),allocatable         :: flux_residue
    real(8),dimension(:,:),allocatable           :: inv_t_xi, t_xi
    real(8),dimension(:,:),allocatable           :: A_pos, A_neg
    real(8),dimension(:,:,:),allocatable         :: A_pos_sys, A_neg_sys
    real(8),dimension(:,:),allocatable           :: inv_teta, teta
    real(8),dimension(:,:),allocatable           :: B_pos, B_neg
    real(8),dimension(:,:,:),allocatable         :: B_pos_sys, B_neg_sys
    real(8),dimension(:,:),allocatable           :: aux_mult
    real(8),dimension(:,:,:),allocatable         :: deltaQ
    real(8),dimension(:,:,:),allocatable         :: E_pos, E_neg, F_pos, F_neg
    real(8),dimension(:,:,:),allocatable         :: main,lower,upper
    real(8),dimension(:,:),allocatable           :: B_sys, aux_deltaQ, deltaQ_til
    integer(4)                                   :: index, index_i, index_j
    real(8)                                      :: eig_pos, eig_neg
    real(8)                                      :: rho_t
    !
allocate(eig(dim))
allocate(flux_residue(imax,jmax,dim))
allocate(deltaQ(imax,jmax,dim))
allocate(aux_mult(dim,dim))
allocate(E_pos(imax,jmax,dim),E_neg(imax,jmax,dim))
allocate(F_pos(imax,jmax,dim),F_neg(imax,jmax,dim))
    ! ******
    ! let is begin using a block tridiagonal algorithm
    ! ******
allocate(main(dim,dim,imax-2),lower(dim,dim,imax-2),upper(dim,dim,imax-2))
allocate(inv_t_xi(dim,dim),t_xi(dim,dim))
allocate(A_pos(dim,dim),A_neg(dim,dim))
allocate(A_pos_sys(dim,dim,imax-2),A_neg_sys(dim,dim,imax-2))
allocate(deltaQ_til(dim,imax-2))
allocate(B_sys(dim,imax-2))
    !
    ! setting the system 
    !
    E_neg = 0.0d0
    F_pos = 0.0d0
    F_neg = 0.0d0
    F_pos = 0.0d0
    call calculate_fluxes(E_pos,E_neg,F_pos,F_neg)    
    !
    !
do j = 2, jmax - 1
    !
    i = 1
    !
    do i = 2, imax - 1
            !
            eig = diag_ksi(U_contravariant(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
            !
            do index = 1, dim 
                call eigen_values_calculate(eig(index),eig_pos,eig_neg)
                diag_pos(index,index) = eig_pos
                diag_neg(index,index) = eig_neg 
            end do 
            !
                inv_t_xi = inv_T_ksi(u(i,j),v(i,j),rho_t,a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
            !
                t_xi     = T_ksi(u(i,j),v(i,j),rho_t,a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
            !
            aux_mult     = matmul(diag_pos,inv_t_xi)
            A_pos        = matmul(t_xi,aux_mult)
            !
            aux_mult     = matmul(diag_neg,inv_t_xi)
            A_neg        = matmul(t_xi,aux_mult)
            !
            do index_i = 1, dim
                do index_j = 1, dim  
                    A_pos_sys(index_i,index_j,i) = A_pos(index_i,index_j)
                    A_neg_sys(index_i,index_j,i) = A_neg(index_i,index_j)
                end do
            end do 
            !
    end do 
    !
    i = imax
    !
    do i = 2, imax - 1
        call residue_flux_vector_splitting_1st(i,j,E_pos,E_neg,F_pos,F_neg,flux_residue)
        B_sys(1,i-1) = -flux_residue(i,j,1) 
        B_sys(2,i-1) = -flux_residue(i,j,2) 
        B_sys(3,i-1) = -flux_residue(i,j,3) 
        B_sys(4,i-1) = -flux_residue(i,j,4) 
    end do
    !
    ! set the matrixes
    !
    do i = 2, imax - 1
        do index_j = 1, dim 
            do index_i = 1, dim 
                main(index_i,index_j,i-1)  = 1.0d0 + A_pos_sys(index_i,index_j,i) - A_neg_sys(index_i,index_j,i) 
                lower(index_i,index_j,i-1) = -A_pos_sys(index_i,index_j,i-1)
                upper(index_i,index_j,i-1) = A_neg_sys(index_i,index_j,i+1)
            end do
        end do 
    end do
    !
    call blktriad(main,lower,upper,dim,imax-2,B_sys,deltaQ_til)
    !
    do i = 2, imax - 1
        deltaQ(i,j,1:dim) = deltaQ_til(1:dim,i-1)
    end do
    !
end do 
!
deallocate(main,lower,upper)
deallocate(inv_t_xi,t_xi)
deallocate(A_pos,A_neg)
deallocate(A_pos_sys,A_neg_sys)
deallocate(deltaQ_til)
deallocate(B_sys)
!
! solving in the eta-direction
!
allocate(aux_deltaQ(dim,jmax-2))
allocate(main(dim,dim,jmax-2),lower(dim,dim,jmax-2),upper(dim,dim,jmax-2))
allocate(inv_teta(dim,dim),teta(dim,dim))
allocate(B_pos(dim,dim),B_neg(dim,dim))
allocate(B_pos_sys(dim,dim,jmax-2),B_neg_sys(dim,dim,jmax-2))
allocate(B_sys(dim,jmax-2))
!
do i = 2, imax - 1
    !
        do j = 2, jmax - 1
        !
        eig = diag_eta(V_contravariant(i,j),a(i,j),eta_x(i,j),eta_y(i,j),dim)
            !
            do index = 1, dim 
                call eigen_values_calculate(eig(index),eig_pos,eig_neg)
                diag_pos(index,index) = eig_pos
                diag_neg(index,index) = eig_neg 
            end do 
            !
                inv_teta = inv_T_eta(u(i,j),v(i,j),rho_t,a(i,j),eta_x(i,j),eta_y(i,j),dim)
            !
                teta     = T_eta(u(i,j),v(i,j),rho_t,a(i,j),eta_x(i,j),eta_y(i,j),dim)
            !
            aux_mult     = matmul(diag_pos,inv_teta)
            B_pos        = matmul(teta,aux_mult)
            !
            aux_mult     = matmul(diag_neg,inv_teta)
            B_neg        = matmul(teta,aux_mult)
            !
            do index_i = 1, dim
                do index_j = 1, dim  
                    B_pos_sys(index_i,index_j,j) = B_pos(index_i,index_j)
                    B_neg_sys(index_i,index_j,j) = B_neg(index_i,index_j)
                end do
            end do 
            !
        end do
        !
    !
    ! set the matrixes
    !
    do j = 2, jmax - 1
        do index_j = 1, dim 
            do index_i = 1, dim 
                main(index_i,index_j,j)  = 1.0d0 + B_pos_sys(index_i,index_j,j) - B_neg_sys(index_i,index_j,j) 
                lower(index_i,index_j,j) = -B_pos_sys(index_i,index_j,j-1)
                upper(index_i,index_j,j) = B_neg_sys(index_i,index_j,j+1)
            end do
        end do 
    end do
    !
    call blktriad(main,lower,upper,dim,jmax-2,B_sys,aux_deltaQ)
    !
end do 
!
deallocate(aux_deltaQ)
deallocate(main,lower,upper)
deallocate(inv_t_xi,t_xi)
deallocate(A_pos,A_neg)
deallocate(A_pos_sys,A_neg_sys)
deallocate(B_sys)
!
!
deallocate(eig)
deallocate(flux_residue)
deallocate(aux_mult)
deallocate(E_pos,E_neg)
deallocate(F_pos,F_neg)
deallocate(deltaQ)
!
!
end subroutine sw_2nd