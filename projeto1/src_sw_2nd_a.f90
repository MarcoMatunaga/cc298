subroutine sw_2nd_a
    use vars
    use vars_sw
    use diagonalization
    use fluxes_pos_neg
    implicit none

    !*******
    ! let is begin using ADI
    !*******
    !*******
    ! create a vector to store the eigenvalues pos and neg
    !*******

    call allocate_vars_sys_ksi
    call allocate_vars_sys_ksi_2nd

    do i = 1, dim
        Identy(i,i,1:imax) = 1.0d0
    end do
   
    ! setting the system 
    
    call residue_flux_vector_splitting_2nd_frontier(E_pos,E_neg,F_pos,F_neg,flux_residue)
    
    do j = 3, jmax - 2
        do i = 3, imax - 2
            !**** colocar a subroutina do residuo com 
            !**** intervalo ao inves de indice por indice
                call residue_flux_vector_splitting_2nd(i,j,E_pos,E_neg,F_pos,F_neg,flux_residue)
            !call residue_flux_vector_splitting_2nd(i,j,E_pos,E_neg,F_pos,F_neg,flux_residue)
            B_sys(1,i-1) = -flux_residue(i,j,1) 
            B_sys(2,i-1) = -flux_residue(i,j,2) 
            B_sys(3,i-1) = -flux_residue(i,j,3) 
            B_sys(4,i-1) = -flux_residue(i,j,4) 
        end do
    end do
    
do j = 2, jmax - 1

    call frontier_xi

    do i = 2, imax - 1

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
            
    end do 
   
    ! set the matrixes
    
    call frontier_xi_2nd
    
    do i = 3, imax - 2
        do index_j = 1, dim 
            do index_i = 1, dim 
                lower_sw(index_i,index_j,i-1)  = 0.50d0*delta_t(i,j)*A_pos_sys(index_i,index_j,i-2)
                lower_sw1(index_i,index_j,i-1) = -2.0d0*delta_t(i,j)*A_pos_sys(index_i,index_j,i-1)
                main_sw(index_i,index_j,i-1)   = Identy(index_i,index_j,i) &
                                                 + 0.50d0*delta_t(i,j)*(3.0d0*A_pos_sys(index_i,index_j,i) &
                                                 -3.0d0*A_neg_sys(index_i,index_j,i))                                                
                upper_sw(index_i,index_j,i-1)  = 2.0d0*delta_t(i,j)*A_neg_sys(index_i,index_j,i+1)
                upper_sw1(index_i,index_j,i-1) = -0.50d0*delta_t(i,j)*A_neg_sys(index_i,index_j,i+2)
            end do
        end do 
    end do

    call penta_block(lower_sw,lower_sw1,main_sw,upper_sw,upper_sw1,B_sys,dim,imax-2,aux_deltaQ_til)

    do i = 2, imax - 1
        deltaQ_til(i,j,1) = aux_deltaQ_til(1,i-1)
        deltaQ_til(i,j,2) = aux_deltaQ_til(2,i-1)
        deltaQ_til(i,j,3) = aux_deltaQ_til(3,i-1)
        deltaQ_til(i,j,4) = aux_deltaQ_til(4,i-1)
    end do
    
end do 

    call deallocate_vars_ksi
    call deallocate_vars_ksi_2nd

! solving in the eta-direction
    
    call allocate_vars_sys_eta
    call allocate_vars_sys_eta_2nd

    do i = 1, dim
        Identy(i,i,1:jmax) = 1.0d0
    end do

do i = 2, imax - 1

        call frontier_eta

        do j = 2, jmax - 1
            
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

        end do

    do j = 2, jmax - 1
        B_sys(1,j-1) = deltaQ_til(i,j,1)
        B_sys(2,j-1) = deltaQ_til(i,j,2)
        B_sys(3,j-1) = deltaQ_til(i,j,3)
        B_sys(4,j-1) = deltaQ_til(i,j,4)
    end do

    ! set the matrixes
    
    call frontier_eta_2nd

    do j = 3, jmax - 2
        do index_j = 1, dim 
            do index_i = 1, dim 
                lower_sw(index_i,index_j,j-1)  = 0.50d0*delta_t(i,j)*B_pos_sys(index_i,index_j,j-2)
                lower_sw1(index_i,index_j,j-1) = -2.0d0*delta_t(i,j)*B_pos_sys(index_i,index_j,j-1)
                main_sw(index_i,index_j,j-1)   = Identy(index_i,index_j,j) &
                                                 + 0.50d0*delta_t(i,j)*(3.0d0*B_pos_sys(index_i,index_j,j) &
                                                 -3.0d0*B_neg_sys(index_i,index_j,j-2))
                upper_sw(index_i,index_j,j-1)  = 2.0d0*delta_t(i,j)*B_neg_sys(index_i,index_j,j+1)
                upper_sw1(index_i,index_j,j-1) = -0.50d0*delta_t(i,j)*B_neg_sys(index_i,index_j,j+2)
            end do
        end do 
    end do

    call penta_block(lower_sw,lower_sw1,main_sw,upper_sw,upper_sw1,B_sys,dim,jmax-2,aux_deltaQ)

    do j = 2, jmax -1
        deltaQ(i,j,1) = aux_deltaQ(1,j-1)
        deltaQ(i,j,2) = aux_deltaQ(2,j-1)
        deltaQ(i,j,3) = aux_deltaQ(3,j-1)
        deltaQ(i,j,4) = aux_deltaQ(4,j-1)
    end do 

    do j = 2, jmax - 1
        Q_barra(i,j,1) = deltaQ(i,j,1) + Q_barra(i,j,1)
        Q_barra(i,j,2) = deltaQ(i,j,2) + Q_barra(i,j,2)
        Q_barra(i,j,3) = deltaQ(i,j,3) + Q_barra(i,j,3)
        Q_barra(i,j,4) = deltaQ(i,j,4) + Q_barra(i,j,4)   
    end do 

end do 
    
        call deallocate_vars_left
        call deallocate_vars_sys_eta_2nd

end subroutine sw_2nd_a