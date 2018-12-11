subroutine  blktriad_pc(maind,lower,upper,id,md,xb,x)

    ! the inverse of a diagonal matrix is one
    ! under each diagonal element

    !|     B(1)    C(1)             | | x(1)  |
    !| A(2)  B(2)    C(2)           | |       |
    !|   A(3)  B(3)                 | |       |
    !|     .     .                  |*|       | = xb[1:mb*3]
    !|       .     .                | |       |
    !|               .       C(mb-1)| |       |
    !|         A(mb)  B(mb)         | |x(n*id)|
    ! id = inner matrices dimension.
    ! md = number matrices.
    ! maind = main diagonal of matrices  format: maind(id,id,md)
    ! lower = lower diagonal of matrices format: lower(id,id,2:md)
    ! upper = upper diagonal of matrices format: maind(id,id,md-1)
    ! xb    = B vector in Ax=B           format: xb(md*id)
    ! x     = x vector in Ax=B           format: x(md*id)
    !
    implicit none
    ! +++ Inputs +++
    !
    ! Scalar input variables.
    integer(kind=4) :: id,md
    ! Main diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: maind
    ! Lower diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: lower
    ! Upper diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: upper
    ! Vector of equalties B in Ax=B.
    real(kind=8), dimension(id,md)    :: xb
    ! Vector of answers x in Ax=B.
    real(kind=8), dimension(id,md)    :: x
    ! ++ Inside variables ++
    !
    ! Scalar variables
    integer(kind=4) :: i,ii
    ! Array of gamma coefficients.
    real(kind=8), dimension(id,id,md) :: gamm
    ! Array of beta coefficients.
    real(kind=8), dimension(id,md) :: beta
    ! Auxiliary arrays.
    real(kind=8), allocatable, dimension(:,:) :: aux_copy
    real(kind=8), allocatable, dimension(:,:) :: aux_mult
    real(kind=8), allocatable, dimension(:,:) :: aux_summ
    real(kind=8), allocatable, dimension(:)   :: aux_dumm
    !
    allocate(aux_mult(id,id))
    allocate(aux_summ(id,id))
    allocate(aux_copy(id,id))
    allocate(aux_dumm(id))
    !--------------------------------------------------------------------------!
    !                  Step 1: BLOCK TRIANGULARIZATION                         !
    !--------------------------------------------------------------------------!
    ! Zero out the auxiliar vectors.
    beta = 0.0d0
    gamm = 0.0d0
    ! Lets first get our first gamma.
    aux_copy = maind(:,:,1)
    !call inv(aux_copy,aux_copy,id)
    do i = 1, id
            aux_copy(i,i) = 1.0d0/aux_copy(i,i)
    end do 
    gamm(:,:,1) = matmul(aux_copy,upper(:,:,1))
    ! Now that we have our first gamma, lets get the rest of then.
    do i = 2, md-1
        aux_mult = 0.0d0
        aux_summ = 0.0d0
        aux_mult = matmul(lower(:,:,i),gamm(:,:,i-1))
            do ii = 1, id
                aux_summ(ii,ii) = maind(ii,ii,i) - aux_mult(ii,ii)
            end do
        !call inv(aux_summ,aux_summ,id)
        do ii = 1, id
            aux_summ(ii,ii) = 1.0d0/aux_summ(ii,ii)
        end do 
        gamm(:,:,i) = matmul(aux_summ,upper(:,:,i))      
    end do
    ! Now that we have our gammas, lets get the betas, starting from the first
    ! ones. Note that now the calls done by the Lapack library will get a bit
    ! more complicated so lets use matmul...
    aux_copy = 0.0d0
    aux_copy = maind(:,:,1)
    !call inv(aux_copy,aux_copy,id)
    do i = 1, id
            aux_copy(i,i) = 1.0d0/aux_copy(i,i)
    end do 
    beta(:,1) = matmul(aux_copy,xb(:,1))
    ! We now have our first beta, lets get the rest.
    do i = 2, md
        aux_mult = 0.0d0
        aux_summ = 0.0d0
        aux_mult(:,:) = matmul(lower(:,:,i),gamm(:,:,i-1))
            do ii = 1, id
                aux_summ(ii,ii) = maind(ii,ii,i) - aux_mult(ii,ii)
            end do
        !call inv(aux_summ,aux_summ,id)
        do ii = 1, id
            aux_summ(ii,ii) = 1.0d0/aux_summ(ii,ii)
        end do 
        aux_dumm(:) = xb(:,i) - matmul(lower(:,:,i),beta(:,i-1))
        beta(:,i) = matmul(aux_summ(:,:),aux_dumm(:))
    end do
    !--------------------------------------------------------------------------!
    !                  Step 2: BACKWARD SWEEP                                  !
    !--------------------------------------------------------------------------!
    ! How cool is that, lets start build our solution vector... iupiiii!
    x = 0.0d0
    x(:,md) = beta(:,md)
    do i = md-1,1,-1
        aux_dumm(:) = matmul(gamm(:,:,i),x(:,i+1))
        do ii = 1, id
            x(ii,i) = beta(ii,i) - aux_dumm(ii)
        end do
    end do
    deallocate(aux_mult)
    deallocate(aux_summ)
    deallocate(aux_copy)
    deallocate(aux_dumm)
end subroutine blktriad_pc
