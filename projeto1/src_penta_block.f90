subroutine penta_block(A_matrix,B_matrix,C_matrix,D_matrix,E_matrix,f_vector,id,nmatrix,x_vector)
    implicit none
    ! the procedure to solve the block pentadiagonals are contained in the following article
    ! An Efficient Implementation of the Thomas-Algorithm for Block Penta-diagonal
    ! Systems on Vector Computers

    !
    integer(4)                                  :: nmatrix, id

    ! input matrices
    real(8), dimension(id,id,nmatrix)           :: C_matrix, D_matrix, E_matrix
    real(8), dimension(id,id,nmatrix)           :: A_matrix
    real(8), dimension(id,id,nmatrix)           :: B_matrix
    real(8), dimension(id,nmatrix)              :: r_vector, x_vector, f_vector
    real(8), dimension(id)                      :: aux_vector_r

    ! 
    real(8), dimension(id,id,-1:nmatrix)        :: Y_matrix, Z_matrix
    real(8), dimension(id,id,1:nmatrix)         :: G_matrix, K_matrix
    real(8), dimension(id,id)                   :: aux_matrix_y, aux_copy

    integer(4)                                  :: ii
    
    !****************
    ! this is made because we are solving in this subroutine two extra points
    ! we change this later
    
    A_matrix(1:id,1:id,1) = 0.0d0
    A_matrix(1:id,1:id,2) = 0.0d0
    B_matrix(1:id,1:id,1) = 0.0d0
    D_matrix(1:id,1:id,nmatrix) = 0.0d0
    E_matrix(1:id,1:id,nmatrix) = 0.0d0
    E_matrix(1:id,1:id,nmatrix-1) = 0.0d0
    
    Y_matrix(1:id,1:id,-1) = 0.0d0
    Y_matrix(1:id,1:id,0) = 0.0d0
    Z_matrix(1:id,1:id,-1) = 0.0d0
    Z_matrix(1:id,1:id,0) = 0.0d0

    do ii = 1, nmatrix
        K_matrix(:,:,ii) = B_matrix(:,:,ii) - matmul(A_matrix(:,:,ii),Y_matrix(:,:,ii-2))
        G_matrix(:,:,ii) = C_matrix(:,:,ii) - matmul(K_matrix(:,:,ii),Y_matrix(:,:,ii-1)) &
                                            - matmul(A_matrix(:,:,ii),Z_matrix(:,:,ii-2))
        ! inv(A,A_inv,m)
        aux_copy = G_matrix(:,:,ii)
        call inv(aux_copy,aux_copy,id)
    
        aux_matrix_y = D_matrix(:,:,ii) - matmul(K_matrix(:,:,ii),Z_matrix(:,:,ii-1))
        aux_vector_r = f_vector(:,ii) - matmul(A_matrix(:,:,ii),r_vector(:,ii-2)) - matmul(K_matrix(:,:,ii),r_vector(:,ii-1))
        Y_matrix(:,:,ii) = matmul(aux_copy,aux_matrix_y)
        Z_matrix(:,:,ii) = matmul(aux_copy,E_matrix(:,:,ii))
        r_vector(:,ii) = matmul(aux_copy,aux_vector_r)
    end do

    ! backward substitution
    x_vector(:,nmatrix)   = r_vector(:,nmatrix)
    x_vector(:,nmatrix-1) = r_vector(:,nmatrix-1) - matmul(Y_matrix(:,:,nmatrix-1),x_vector(:,nmatrix))
    do ii = nmatrix-2,1,-1
        x_vector(:,ii) = r_vector(:,ii) - matmul(Y_matrix(:,:,ii),x_vector(:,ii+1)) - matmul(Z_matrix(:,:,ii),x_vector(:,ii+2))
    end do

end subroutine penta_block