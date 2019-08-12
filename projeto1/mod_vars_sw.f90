module vars_sw
    use vars
    use diagonalization
    implicit none
    real(8),dimension(:),allocatable             :: eig
    real(8),dimension(:,:),allocatable           :: diag_pos, diag_neg
    real(8),dimension(:,:),allocatable           :: inv_t_xi, t_xi
    real(8),dimension(:,:),allocatable           :: A_pos, A_neg
    real(8),dimension(:,:),allocatable           :: inv_teta, teta
    real(8),dimension(:,:),allocatable           :: B_pos, B_neg
    real(8),dimension(:,:),allocatable           :: aux_mult
    real(8),dimension(:,:),allocatable           :: B_sys, aux_deltaQ, aux_deltaQ_til
    real(8),dimension(:,:),allocatable           :: u,v,rho,p,a
    real(8),dimension(:,:,:),allocatable         :: flux_residue, Identy
    real(8),dimension(:,:,:),allocatable         :: A_pos_sys, A_neg_sys
    real(8),dimension(:,:,:),allocatable         :: B_pos_sys, B_neg_sys
    real(8),dimension(:,:,:),allocatable         :: deltaQ, deltaQ_til
    real(8),dimension(:,:,:),allocatable         :: main_sw,lower_sw,upper_sw
    real(8),dimension(:,:,:),allocatable         :: lower_sw1,upper_sw1
    integer(4)                                   :: index, index_i, index_j
    real(8)                                      :: eig_pos, eig_neg
contains
    
    subroutine allocate_vars_sw
        implicit none
        allocate(eig(dim))
        allocate(diag_pos(dim,dim),diag_neg(dim,dim))
        allocate(flux_residue(imax,jmax,dim))
        allocate(deltaQ(imax,jmax,dim))
        allocate(deltaQ_til(imax,jmax,dim))
        allocate(aux_mult(dim,dim))
        allocate(u(imax,jmax),v(imax,jmax),rho(imax,jmax),p(imax,jmax),a(imax,jmax))
    end subroutine allocate_vars_sw

    subroutine allocate_vars_sys_ksi
        implicit none
        allocate(main_sw(dim,dim,imax-2),lower_sw(dim,dim,imax-2),upper_sw(dim,dim,imax-2))
        allocate(inv_t_xi(dim,dim),t_xi(dim,dim))
        allocate(A_pos(dim,dim),A_neg(dim,dim))
        allocate(A_pos_sys(dim,dim,imax),A_neg_sys(dim,dim,imax),Identy(dim,dim,imax))
        allocate(aux_deltaQ_til(dim,imax-2))
        allocate(B_sys(dim,imax-2))
        diag_neg       = 0.0d0
        diag_pos       = 0.0d0
        lower_sw       = 0.0d0
        main_sw        = 0.0d0 
        upper_sw       = 0.0d0
        B_sys          = 0.0d0
        aux_deltaQ_til = 0.0d0
        deltaQ_til     = 0.0d0
        deltaQ         = 0.0d0  
        A_pos_sys      = 0.0d0
        A_pos          = 0.0d0
        A_neg_sys      = 0.0d0
        A_neg          = 0.0d0
        Identy         = 0.0d0
        inv_t_xi       = 0.0d0
        t_xi           = 0.0d0
        flux_residue   = 0.0d0
        aux_mult       = 0.0d0
    end subroutine allocate_vars_sys_ksi

    subroutine allocate_vars_sys_ksi_2nd
        implicit none
        allocate(lower_sw1(dim,dim,imax-2),upper_sw1(dim,dim,imax-2))
        lower_sw1     = 0.0d0
        upper_sw1     = 0.0d0
    end subroutine allocate_vars_sys_ksi_2nd

    subroutine deallocate_vars_ksi
        implicit none
        deallocate(main_sw,lower_sw,upper_sw)
        deallocate(inv_t_xi,t_xi)
        deallocate(A_pos,A_neg)
        deallocate(A_pos_sys,A_neg_sys)
        deallocate(aux_deltaQ_til)
        deallocate(B_sys)
        deallocate(Identy)
    end subroutine deallocate_vars_ksi

    subroutine deallocate_vars_ksi_2nd
        implicit none
        deallocate(lower_sw1,upper_sw1)
    end subroutine deallocate_vars_ksi_2nd

    subroutine allocate_vars_sys_eta
        implicit none
        allocate(aux_deltaQ(dim,jmax-2))
        allocate(main_sw(dim,dim,jmax-2),lower_sw(dim,dim,jmax-2),upper_sw(dim,dim,jmax-2))
        allocate(inv_teta(dim,dim),teta(dim,dim))
        allocate(B_pos(dim,dim),B_neg(dim,dim),Identy(dim,dim,jmax))
        allocate(B_pos_sys(dim,dim,jmax),B_neg_sys(dim,dim,jmax))
        allocate(B_sys(dim,jmax-2))

        Identy     = 0.0d0
        aux_deltaQ = 0.0d0
        lower_sw   = 0.0d0
        main_sw    = 0.0d0
        upper_sw   = 0.0d0
        B_sys      = 0.0d0
        B_pos_sys  = 0.0d0
        B_neg_sys  = 0.0d0  
        teta       = 0.0d0
        inv_teta   = 0.0d0 
    end subroutine allocate_vars_sys_eta

    subroutine allocate_vars_sys_eta_2nd
        implicit none
        allocate(lower_sw1(dim,dim,jmax-2),upper_sw1(dim,dim,jmax-2))
        lower_sw1   = 0.0d0
        upper_sw1   = 0.0d0
    end subroutine allocate_vars_sys_eta_2nd

    subroutine deallocate_vars_sys_eta_2nd
        implicit none
        deallocate(lower_sw1,upper_sw1)
    end subroutine deallocate_vars_sys_eta_2nd
    
    subroutine deallocate_vars_left
        implicit none
        deallocate(aux_deltaQ)
        deallocate(main_sw,lower_sw,upper_sw)
        deallocate(B_sys,B_pos,B_neg)
        deallocate(B_pos_sys,B_neg_sys)
        deallocate(Identy)
        
        deallocate(inv_teta,teta)
        deallocate(eig)
        deallocate(diag_pos,diag_neg)
        deallocate(flux_residue)
        deallocate(aux_mult)
        deallocate(deltaQ)
        deallocate(deltaQ_til)
        deallocate(u,v,rho,p,a)        
    end subroutine deallocate_vars_left

    subroutine frontier_xi
        implicit none
        
        i = 1 

        eig = diag_ksi(U_contravariant(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),dim)

            do index = 1, dim 
                eig_pos = 0.50d0*(eig(index) + abs(eig(index)))
                eig_neg = 0.50d0*(eig(index) - abs(eig(index)))
                diag_pos(index,index) = eig_pos
                diag_neg(index,index) = eig_neg 
            end do 
            
            inv_t_xi = inv_T_ksi(u(i,j),v(i,j),rho(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
            
            t_xi     = T_ksi(u(i,j),v(i,j),rho(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
            
            aux_mult     = matmul(diag_pos,inv_t_xi)
            A_pos        = matmul(t_xi,aux_mult)
            
            aux_mult     = matmul(diag_neg,inv_t_xi)
            A_neg        = matmul(t_xi,aux_mult)
            
            do index_j = 1, dim
                do index_i = 1, dim  
                    A_pos_sys(index_i,index_j,i) = A_pos(index_i,index_j)
                    A_neg_sys(index_i,index_j,i) = A_neg(index_i,index_j)
                end do
            end do 

        i = imax 

        eig = diag_ksi(U_contravariant(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),dim)

            do index = 1, dim 
                eig_pos = 0.50d0*(eig(index) + abs(eig(index)))
                eig_neg = 0.50d0*(eig(index) - abs(eig(index)))
                diag_pos(index,index) = eig_pos
                diag_neg(index,index) = eig_neg 
            end do 
            
            inv_t_xi = inv_T_ksi(u(i,j),v(i,j),rho(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
            
            t_xi     = T_ksi(u(i,j),v(i,j),rho(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),dim)
            
            aux_mult     = matmul(diag_pos,inv_t_xi)
            A_pos        = matmul(t_xi,aux_mult)
            
            aux_mult     = matmul(diag_neg,inv_t_xi)
            A_neg        = matmul(t_xi,aux_mult)
            
            do index_j = 1, dim
                do index_i = 1, dim  
                    A_pos_sys(index_i,index_j,i) = A_pos(index_i,index_j)
                    A_neg_sys(index_i,index_j,i) = A_neg(index_i,index_j)
                end do
            end do 

    end subroutine frontier_xi

    subroutine frontier_eta
        implicit none
        
        j = 1 

        eig = diag_eta(V_contravariant(i,j),a(i,j),eta_x(i,j),eta_y(i,j),dim)

            do index = 1, dim 
                eig_pos = 0.50d0*(eig(index) + abs(eig(index)))
                eig_neg = 0.50d0*(eig(index) - abs(eig(index)))
                diag_pos(index,index) = eig_pos
                diag_neg(index,index) = eig_neg 
            end do 

            inv_teta = inv_T_eta(u(i,j),v(i,j),rho(i,j),a(i,j),eta_x(i,j),eta_y(i,j),dim)

            teta     = T_eta(u(i,j),v(i,j),rho(i,j),a(i,j),eta_x(i,j),eta_y(i,j),dim)

            aux_mult     = matmul(diag_pos,inv_teta)
            B_pos        = matmul(teta,aux_mult)

            aux_mult     = matmul(diag_neg,inv_teta)
            B_neg        = matmul(teta,aux_mult)

            do index_i = 1, dim
                do index_j = 1, dim  
                    B_pos_sys(index_i,index_j,j) = B_pos(index_i,index_j)
                    B_neg_sys(index_i,index_j,j) = B_neg(index_i,index_j)
                end do
            end do 
        
        j = jmax 

            do index = 1, dim 
                eig_pos = 0.50d0*(eig(index) + abs(eig(index)))
                eig_neg = 0.50d0*(eig(index) - abs(eig(index)))
                diag_pos(index,index) = eig_pos
                diag_neg(index,index) = eig_neg 
            end do 

            inv_teta = inv_T_eta(u(i,j),v(i,j),rho(i,j),a(i,j),eta_x(i,j),eta_y(i,j),dim)

            teta     = T_eta(u(i,j),v(i,j),rho(i,j),a(i,j),eta_x(i,j),eta_y(i,j),dim)

            aux_mult     = matmul(diag_pos,inv_teta)
            B_pos        = matmul(teta,aux_mult)

            aux_mult     = matmul(diag_neg,inv_teta)
            B_neg        = matmul(teta,aux_mult)

            do index_i = 1, dim
                do index_j = 1, dim  
                    B_pos_sys(index_i,index_j,j) = B_pos(index_i,index_j)
                    B_neg_sys(index_i,index_j,j) = B_neg(index_i,index_j)
                end do
            end do 

    end subroutine frontier_eta

    subroutine frontier_xi_2nd
        implicit none
        
        do index_j = 1, dim 
            do index_i = 1, dim 
        
            i = 2
                lower_sw(index_i,index_j,i-1)  = 0.0d0
                lower_sw1(index_i,index_j,i-1) = -2.0d0*delta_t(i,j)*A_pos_sys(index_i,index_j,i-1)
                main_sw(index_i,index_j,i-1)   = Identy(index_i,index_j,i) & 
                                                 + 0.50d0*delta_t(i,j)*(3.0d0*A_pos_sys(index_i,index_j,i) &
                                                 -3.0d0*A_neg_sys(index_i,index_j,i))                                                
                upper_sw(index_i,index_j,i-1)  = 2.0d0*delta_t(i,j)*A_neg_sys(index_i,index_j,i+1)
                upper_sw1(index_i,index_j,i-1) = -0.50d0*delta_t(i,j)*A_neg_sys(index_i,index_j,i+2)
    
            i = imax - 1
                lower_sw(index_i,index_j,i-1)  = 0.50d0*delta_t(i,j)*A_pos_sys(index_i,index_j,i-2)
                lower_sw1(index_i,index_j,i-1) = -2.0d0*delta_t(i,j)*A_pos_sys(index_i,index_j,i-1)
                main_sw(index_i,index_j,i-1)   = Identy(index_i,index_j,i) & 
                                                 + 0.50d0*delta_t(i,j)*(3.0d0*A_pos_sys(index_i,index_j,i) &
                                                 -3.0d0*A_neg_sys(index_i,index_j,i))                                                
                upper_sw(index_i,index_j,i-1)  = 2.0d0*delta_t(i,j)*A_neg_sys(index_i,index_j,i+1)
                upper_sw1(index_i,index_j,i-1) = 0.0d0
            end do 
        end do 

    end subroutine frontier_xi_2nd

    subroutine frontier_eta_2nd
        implicit none
        
        do index_j = 1, dim 
            do index_i = 1, dim 

            j = 2
                lower_sw(index_i,index_j,j-1)  = 0.0d0
                lower_sw1(index_i,index_j,j-1) = -2.0d0*delta_t(i,j)*B_pos_sys(index_i,index_j,j-1)
                main_sw(index_i,index_j,j-1)   = Identy(index_i,index_j,j) & 
                                                 + 0.50d0*delta_t(i,j)*(3.0d0*B_pos_sys(index_i,index_j,j) &
                                                 -3.0d0*B_neg_sys(index_i,index_j,j))
                upper_sw(index_i,index_j,j-1)  = 2.0d0*delta_t(i,j)*B_neg_sys(index_i,index_j,j+1)
                upper_sw1(index_i,index_j,j-1) = -0.50d0*delta_t(i,j)*B_neg_sys(index_i,index_j,j+2)

            j = jmax - 1
                lower_sw(index_i,index_j,j-1)  = 0.50d0*delta_t(i,j)*B_pos_sys(index_i,index_j,j-2)
                lower_sw1(index_i,index_j,j-1) = -2.0d0*delta_t(i,j)*B_pos_sys(index_i,index_j,j-1)
                main_sw(index_i,index_j,j-1)   = Identy(index_i,index_j,j) & 
                                                 + 0.50d0*delta_t(i,j)*(3.0d0*B_pos_sys(index_i,index_j,j) &
                                                 -3.0d0*B_neg_sys(index_i,index_j,j))
                upper_sw(index_i,index_j,j-1)  = 2.0d0*delta_t(i,j)*B_neg_sys(index_i,index_j,j+1)
                upper_sw1(index_i,index_j,j-1) = 0.0d0
            end do 
        end do 

    end subroutine frontier_eta_2nd
end module vars_sw