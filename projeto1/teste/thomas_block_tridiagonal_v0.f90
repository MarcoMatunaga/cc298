subroutine thomas_block_tridiagonal(x_dir,y_dir,size,lower,main,upper,sol,res)
    use vars
    implicit none
!
! the system that will be solved is denoted by A X = B 
! where A is the block tridiagonal, x the vector of 'solutions', and B the vector
! of 'results'
!
! the lower diagonal is denote by a_matrix, 
! the main diagonal is denoted by b_matrix,
! the upper diagonal is denoted by c_matrix.
!
! lower, main and upper are matrices
!
! sol and res are matrices as well just to make
! easier the implementation (see pulliam to understand)
!
! size is the number of elemensts in the main diagonal
!
integer size, i_sub, j_sub, diag, num_call
integer x_dir, y_dir, index_dir
double precision sol(imax,jmax,dim), res(imax,jmax,dim)
double precision lower(size,dim,dim), main(size,dim,dim), upper(size,dim,dim)
!
! declare the variables used to solve the block tridiagonal system
! this algorithm follows the pdf in folder project1 (section4.pdf)
!
double precision b_matrix(dim,dim), inv_matrix(dim,dim), c_matrix(dim,dim), a_matrix(dim,dim)
double precision gama_matrix(size,dim,dim), aux(dim,dim), aux1(dim,dim), aux2(dim,dim), aux3(dim,dim)
!
! vectors of the linear system
!
double precision beta(size,dim), x(size,dim), d(size,dim)
double precision beta_aux(dim), x_aux(dim), d_aux(dim), aux_v(dim), aux_v1(dim)
!
! let is start
!
b_matrix = 0.0d0
c_matrix = 0.0d0
inv_matrix = 0.0d0
a_matrix = 0.0d0
gama_matrix = 0.0d0
aux = 0.0d0
aux1 = 0.0d0
aux2 = 0.0d0
aux3 = 0.0d0
aux_v = 0.0d0
aux_v1 = 0.0d0
beta_aux = 0.0d0
d_aux = 0.0d0
x_aux = 0.0d0
!
! gama(1) = (inverse of b1) * (c1)
!
do j_sub = 1, dim
    do i_sub = 1, dim
        b_matrix(i_sub,j_sub) = main(1,i_sub,j_sub)
        c_matrix(i_sub,j_sub) = upper(1,i_sub,j_sub)
    end do
end do
call inverse(b_matrix,inv_matrix,dim)
aux = matmul(inv_matrix,c_matrix)  
do j_sub = 1, dim
    do i_sub = 1, dim
        gama_matrix(1,i_sub,j_sub) = aux(i_sub,j_sub)
    end do 
end do
!
!
aux = 0.0d0
inv_matrix = 0.0d0
!
! solve from the second diagonal to the size - 1
! diag = 2, size - 1
!
do diag = 2, size - 1
    !
    !
    do j_sub = 1, dim
        do i_sub = 1, dim
            aux1(i_sub,j_sub) = lower(diag,i_sub,j_sub)
            aux2(i_sub,j_sub) = gama_matrix(diag-1,i_sub,j_sub)
        end do 
    end do
    !
    !
    aux3 = matmul(aux1,aux2)
    !
    !
    aux1 = 0.0d0
    !
    !
    do j_sub = 1, dim
        do i_sub = 1, dim
            aux(i_sub,j_sub) = main(diag,i_sub,j_sub) - aux3(i_sub,j_sub)
        end do 
    end do
    !
    !
    aux3 = 0.0d0
    !
    !
    call inverse(aux,inv_matrix,dim)
    !
    !
    do j_sub = 1, dim
        do i_sub = 1, dim
            c_matrix(i_sub,j_sub) = upper(diag,i_sub,j_sub)
        end do
    end do
    !
    !
    aux1 = matmul(inv_matrix,c_matrix)
    !
    !
    do j_sub = 1, dim
        do i_sub = 1, dim
            gama_matrix(diag,i_sub,j_sub) = aux1(i_sub,j_sub)
        end do 
    end do
    !
    !
end do
!
!
aux = 0.0d0
aux1 = 0.0d0
inv_matrix = 0.0d0
!
!
diag = 1
!
!
do j_sub = 1, dim
    do i_sub = 1, dim
        b_matrix(i_sub,j_sub) = main(diag,i_sub,j_sub)
    end do
end do
!
call inverse(b_matrix,inv_matrix,dim)
!**************************************
! too many ifs one by subroutine call
!**************************************
if (x_dir == 1) then
    do i_sub = 1, dim 
        d_aux(i_sub) = res(diag,y_dir,i_sub) 
    end do
else
    do j_sub = 1, dim 
        d_aux(j_sub) = res(x_dir,diag,j_sub) 
    end do
end if
!
!
beta_aux = matmul(inv_matrix,d_aux)
!
!
do i_sub = 1, dim
    beta(diag,i_sub) = beta_aux(i_sub)
end do
!
!
beta_aux = 0.0d0
!
!
do diag = 2, size
    !
    !
    do j_sub = 1, dim
        do i_sub = 1, dim
            b_matrix(i_sub,j_sub) = main(diag,i_sub,j_sub) 
            aux(i_sub,j_sub) = lower(diag,i_sub,j_sub)
            aux1(i_sub,j_sub) = gama_matrix(diag-1,i_sub,j_sub)
        end do
    end do
    !
    !
    aux2 = matmul(aux,aux1)
    aux = 0.0d0
    aux1 = 0.0d0
    !
    !
    do j_sub = 1, dim
        do i_sub = 1, dim
            aux(i_sub,j_sub) = b_matrix(i_sub,j_sub) - aux2(i_sub,j_sub)
        end do 
    end do
    !
    !
    inv_matrix = 0.0d0
    call inverse(aux,inv_matrix,dim)
    aux = 0.0d0
    aux2 = 0.0d0
    !
    !
    if ( x_dir == 1 ) then
        do i_sub = 1, dim
            d_aux(i_sub) = res(diag,y_dir,i_sub)
        end do
    else
        do j_sub = 1, dim
            d_aux(j_sub) = res(x_dir,diag,j_sub)
        end do
    end if
    !
    !
    do j_sub = 1, dim
        do i_sub = 1, dim
            a_matrix(i_sub,j_sub) = lower(diag,i_sub,j_sub)
        end do
    end do
    !
    !
    do i_sub = 1, dim
        beta_aux(i_sub) = beta(diag - 1,i_sub)
    end do
    !
    ! ****************
    ! be careful with this operation
    ! ****************
    !
    aux_v = d_aux - matmul(a_matrix,beta_aux)
    aux_v1 = matmul(inv_matrix,aux_v)
    !
    !
    do i_sub = 1, dim
        beta(diag,i_sub) = aux_v1(i_sub)
    end do
    !
    !
    aux_v = 0.0d0
    aux_v1 = 0.0d0
    beta_aux = 0.0d0
    !
    !
end do
!
! back substitution
!
diag = size
!
!
do i_sub = 1, dim
    x(diag,i_sub) = beta(diag,i_sub)
end do
!
!
do diag = size - 1, 1,-1
    !
    !
    do i_sub = 1, dim
        beta_aux(i_sub) = beta(diag,i_sub) 
        x_aux(i_sub) = x(diag+1,i_sub)
    end do
    !
    !
    do j_sub = 1, dim
        do i_sub = 1, dim
            aux(i_sub,j_sub) = gama_matrix(diag,i_sub,j_sub)
        end do 
    end do
    !
    ! ****************
    ! be careful with this operation
    ! ****************
    !
    aux_v = beta_aux - matmul(aux,x_aux)
    !
    !
    do i_sub = 1, dim
        x(diag,i_sub) = aux_v(i_sub) 
    end do
    !
    !
    aux_v = 0.d0
    aux = 0.0d0
    !
    !
end do
!
! i_sub and j_sub change their usual value only here
! because we gonna update the value of the solution matrix
!
if ( x_dir == 1 ) then 
    !
    !
    do i_sub = 1, imax
        sol(i_sub,y_dir,1) = x(i_sub,1)
        sol(i_sub,y_dir,2) = x(i_sub,2)
        sol(i_sub,y_dir,3) = x(i_sub,3)
        sol(i_sub,y_dir,4) = x(i_sub,4)
    end do 
    !
    !
else
    !
    !
    do j_sub = 1, jmax
        sol(x_dir,j_sub,1) = x(j_sub,1)
        sol(x_dir,j_sub,2) = x(j_sub,2)
        sol(x_dir,j_sub,3) = x(j_sub,3)
        sol(x_dir,j_sub,4) = x(j_sub,4)
    end do
    !
    !
end if
!
!
!
end subroutine thomas_block_tridiagonal
