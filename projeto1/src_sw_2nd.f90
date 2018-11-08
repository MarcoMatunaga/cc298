subroutine sw_2nd
    use vars
    use vars_sw
    use diagonalization
    use fluxes_pos_neg
    implicit none

    ! ******
    ! let is begin using ADI
    ! ******

    call allocate_vars_sw
    
    call allocate_vars_sys_ksi
    do i = 1, dim
        Identy(i,i,1:imax) = 1.0d0
    end do

    call calculate_fluxes(E_pos,E_neg,F_pos,F_neg)
    
    ! setting the system 
    
do j = 2, jmax - 1
    
    do i = 2, imax - 1
            eig = diag_ksi(U_contravariant(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),dim)

            do index = 1, dim 
                call eigen_values_calculate(eig(index),eig_pos,eig_neg)
                diag_pos(index,index) = eig_pos
                diag_neg(index,index) = eig_neg 
            end do 
            !
            rho_t = Q_barra(i,j,1)/metric_jacobian(i,j)
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
            do index_j = 1, dim
                do index_i = 1, dim  
                    A_pos_sys(index_i,index_j,i) = delta_t(i,j)*A_pos(index_i,index_j)
                    A_neg_sys(index_i,index_j,i) = delta_t(i,j)*A_neg(index_i,index_j)
                end do
            end do 
    end do 
    
    do i = 2, imax - 1
        !**** colocar a subroutina do residuo com 
        !**** intervalo ao inves de indice por indice
        call residue_flux_vector_splitting_1st(i,j,E_pos,E_neg,F_pos,F_neg,flux_residue)
        B_sys(1,i-1) = -flux_residue(i,j,1) 
        B_sys(2,i-1) = -flux_residue(i,j,2) 
        B_sys(3,i-1) = -flux_residue(i,j,3) 
        B_sys(4,i-1) = -flux_residue(i,j,4) 
    end do
    
    ! set the matrixes

    do i = 2, imax - 1
        do index_j = 1, dim 
            do index_i = 1, dim 
                main_sw(index_i,index_j,i-1)  = Identy(index_i,index_j,i) + A_pos_sys(index_i,index_j,i) &
                                                -A_neg_sys(index_i,index_j,i) 
                lower_sw(index_i,index_j,i-1) = -A_pos_sys(index_i,index_j,i-1)
                upper_sw(index_i,index_j,i-1) = A_neg_sys(index_i,index_j,i+1)
            end do
        end do 
    end do

    call blktriad(main_sw,lower_sw,upper_sw,dim,imax-2,B_sys,aux_deltaQ_til)

    do i = 2, imax - 1
        deltaQ_til(i,j,1) = aux_deltaQ_til(1,i-1)
        deltaQ_til(i,j,2) = aux_deltaQ_til(2,i-1)
        deltaQ_til(i,j,3) = aux_deltaQ_til(3,i-1)
        deltaQ_til(i,j,4) = aux_deltaQ_til(4,i-1)
    end do
    
end do 

    call deallocate_vars_ksi

! solving in the eta-direction
    
    call allocate_vars_sys_eta

    do i = 1, dim
        Identy(i,i,1:jmax) = 1.0d0
    end do

do i = 2, imax - 1

        do j = 2, jmax - 1
            eig = diag_eta(V_contravariant(i,j),a(i,j),eta_x(i,j),eta_y(i,j),dim)

            do index = 1, dim 
                call eigen_values_calculate(eig(index),eig_pos,eig_neg)
                diag_pos(index,index) = eig_pos
                diag_neg(index,index) = eig_neg 
            end do 

            rho_t = Q_barra(i,j,1)/metric_jacobian(i,j)
            inv_teta = inv_T_eta(u(i,j),v(i,j),rho_t,a(i,j),eta_x(i,j),eta_y(i,j),dim)

            teta     = T_eta(u(i,j),v(i,j),rho_t,a(i,j),eta_x(i,j),eta_y(i,j),dim)

            aux_mult     = matmul(diag_pos,inv_teta)
            B_pos        = matmul(teta,aux_mult)

            aux_mult     = matmul(diag_neg,inv_teta)
            B_neg        = matmul(teta,aux_mult)

            do index_i = 1, dim
                do index_j = 1, dim  
                    B_pos_sys(index_i,index_j,j) = delta_t(i,j)*B_pos(index_i,index_j)
                    B_neg_sys(index_i,index_j,j) = delta_t(i,j)*B_neg(index_i,index_j)
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

    do j = 2, jmax - 1
        do index_j = 1, dim 
            do index_i = 1, dim 
                main_sw(index_i,index_j,j-1)  = Identy(index_i,index_j,j) + B_pos_sys(index_i,index_j,j) &
                                                -B_neg_sys(index_i,index_j,j) 
                lower_sw(index_i,index_j,j-1) = -B_pos_sys(index_i,index_j,j-1)
                upper_sw(index_i,index_j,j-1) = B_neg_sys(index_i,index_j,j+1)
            end do
        end do 
    end do

    call blktriad(main_sw,lower_sw,upper_sw,dim,jmax-2,B_sys,aux_deltaQ)

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

        ! if (iter == 8) then 
        !     open(997,file='sw')
        !     do j = 2, jmax - 1
        !         do i = 2, imax - 1 
        !             write(997,*) 'pos1',i,j,Q_barra(i,j,1)
        !             write(997,*) 'pos2',i,j,Q_barra(i,j,2)
        !             write(997,*) 'pos3',i,j,Q_barra(i,j,3)
        !             write(997,*) 'pos4',i,j,Q_barra(i,j,4)
        !         end do
        !     end do 
        !     close(997)    
        ! end if

        ! if (iter == 0) then 
        !     open (996,file='t_matrices_sw') 
        !     do j = 2, jmax - 1
        !         do i = 2, imax - 1 
        !             write(996,*)  'ksi',U_contravariant(i,j),a(i,j),ksi_x(i,j),ksi_y(i,j),&
        !                                 u(i,j),v(i,j),Q_barra(i,j,1)/metric_jacobian(i,j)
        !             write(996,*)  'eta',V_contravariant(i,j),a(i,j),eta_x(i,j),eta_y(i,j),&
        !                                 u(i,j),v(i,j),Q_barra(i,j,1)/metric_jacobian(i,j)
        !         end do
        !     end do 
        !     close(996) 
        ! end if

        call deallocate_vars_left
        
end subroutine sw_2nd

