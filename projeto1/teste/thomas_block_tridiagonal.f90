program testando

    ! This program main goal is to test the block tridiagonal solver developed
    ! for CFD applications. There are a lot of good optimizations that can be
    ! done here but this subroutine works perfectly.

    ! Compile the code with gfortran [this file] -lblas -llapack
    ! Run ./a.out
    ! The screen output should be the correct answer hard coded below.

    implicit none

    integer(4) :: i,j, jmax, size_matrix
    real(8), dimension(2,2,3) :: A
    real(8), dimension(2,2,3) :: B,C
    real(8), dimension(2,3) :: b_sys
    real(8), dimension(2,3) :: x

    size_matrix = 2
    !
    ! Define the matrix to be tested. It's a block matrix.
    !

    do i = 1,3
        A(1,1,i) = 4.0d0
        A(1,2,i) = 3.0d0
        A(2,1,i) = 3.0d0
        A(2,2,i) = 2.0d0
    end do

    do i = 1,3

        B(1,1,i) = 8.0d0
        B(1,2,i) = 6.0d0
        B(2,1,i) = 6.0d0
        B(2,2,i) = 4.0d0

        C(1,1,i) = 16.0d0
        C(1,2,i) = 9.0d0
        C(2,1,i) = 9.0d0
        C(2,2,i) = 6.0d0
    end do

    !
    ! Define the B vector in AX = B system.
    !

    do i = 1, 2
        do j = 1, 3
            b_sys(i,j) = 3.0d0
        end do
    end do

    !
    ! Call solver.
    !

    call thomas_block_tridiagonal(A,B,C,size_matrix,3,b_sys,x)


    !
    ! Hard coded answer.
    !

    ! write(*,*) "x(1) = 5.8571429"
    ! write(*,*) "x(2) = -8.9220779"
    ! write(*,*) "x(3) = 0.5714286"
    ! write(*,*) "x(4) = -0.3116883"
    ! write(*,*) "x(5) = 1.8571429"
    ! write(*,*) "x(6) = -2.3766234"

    ! write(*,*) "Calculated..." 

    ! do j = 1, 3
    !     do i = 1, 2
    !         write(*,*) "x(",i,") = ",x(i,j)
    !     end do
    ! end do

end program testando
!
!
!
subroutine thomas_block_tridiagonal(lower,main,upper,dim,diag,res,sol)
    implicit none
    !
    ! input values and output sol
    !
    integer dim, diag
    double precision lower(dim,dim,diag), main(dim,dim,diag), upper(dim,dim,diag)
    double precision sol(dim,diag), res(dim,diag)
    !
    !
    double precision gama_matrix(dim,dim,diag)
    double precision beta(dim,diag)
    !
    ! auxiliary vectors 
    !
    double precision a_matrix(dim,dim), b_matrix(dim,dim), c_matrix(dim,dim)
    double precision x_v(dim), d_v(dim)
    !
    !
    double precision gama_aux(dim,dim), inv_matrix(dim,dim)
    double precision beta_aux(dim)
    double precision aux(dim,dim), aux_v(dim)
    !
    !
    integer i_sub, j_sub, i_diag
    !
    !
    i_diag = 1
    !
    !
    do j_sub = 1, dim
        do i_sub = 1, dim
            b_matrix(i_sub,j_sub) = main(i_sub,j_sub,i_diag)
            c_matrix(i_sub,j_sub) = upper(i_sub,j_sub,i_diag)
        end do 
    end do 
    call inverse(b_matrix,inv_matrix,dim)
        ! do j_sub = 1, dim
        !     do i_sub = 1, dim
        !         write(*,*) inv_matrix(i_sub,j_sub)
        !     end do 
        ! end do
    !
    !
    aux = matmul(inv_matrix,c_matrix)
        ! do j_sub = 1, dim
        !     do i_sub = 1, dim
        !         write(*,*) aux(i_sub,j_sub)
        !     end do 
        ! end do
    do j_sub = 1, dim
        do i_sub = 1, dim
            gama_matrix(i_sub,j_sub,i_diag) = aux(i_sub,j_sub)
        end do 
    end do
    aux = 0.0d0
    inv_matrix = 0.0d0
    ! do i_diag = 1, diag
    !     do j_sub = 1, dim
    !         do i_sub = 1, dim
    !             write(*,*) gama_matrix(i_sub,j_sub,i_diag)
    !         end do 
    !     end do
    ! end do
    !
    !
    !
    do i_diag = 2, diag - 1
        !
        !
        do j_sub = 1, dim
            do i_sub = 1, dim
                gama_aux(i_sub,j_sub) = gama_matrix(i_sub,j_sub,i_diag-1)
                a_matrix(i_sub,j_sub) = lower(i_sub,j_sub,i_diag)
                b_matrix(i_sub,j_sub) = main(i_sub,j_sub,i_diag)           
                c_matrix(i_sub,j_sub) = upper(i_sub,j_sub,i_diag)
            end do 
        end do
        !
        !
        aux = b_matrix - matmul(a_matrix,gama_aux)
        do j_sub = 1, dim
            do i_sub = 1, dim
                write(*,*) i_sub,j_sub,aux(i_sub,j_sub)
            end do 
        end do
        call inverse(aux,inv_matrix,dim)
        ! do j_sub = 1, dim
        !     do i_sub = 1, dim
        !         write(*,*) inv_matrix(i_sub,j_sub)
        !     end do 
        ! end do
        aux = 0.0d0
        aux = matmul(inv_matrix,c_matrix)
        ! do j_sub = 1, dim
        !     do i_sub = 1, dim
        !         write(*,*) aux(i_sub,j_sub)
        !     end do 
        ! end do
        do j_sub = 1, dim
            do i_sub = 1, dim
                gama_matrix(i_sub,j_sub,i_diag) = aux(i_sub,j_sub)
            end do 
        end do 
        !
        !
        inv_matrix = 0.0d0
    end do
    !
    !
    i_diag = 1
    !
    !
    do j_sub = 1, dim
        do i_sub = 1, dim 
            b_matrix(i_sub,j_sub) = main(i_sub,j_sub,i_diag)
        end do 
    end do
    !
    !
    call inverse(b_matrix,inv_matrix,dim)
    do i_sub = 1, dim 
        d_v(i_sub) = res(i_sub,i_diag) 
    end do 
    !
    !
    aux_v = matmul(inv_matrix,d_v)
    do i_sub = 1, dim 
        beta(i_sub,i_diag) = aux_v(i_sub)
    end do 
    !
    !
    inv_matrix = 0.0d0
    !
    !
    do i_diag = 2, diag
        !
        !
        do j_sub = 1, dim 
            do i_sub = 1, dim 
                b_matrix(i_sub,j_sub) = main(i_sub,j_sub,i_diag) 
                a_matrix(i_sub,j_sub) = lower(i_sub,j_sub,i_diag)
                gama_aux(i_sub,j_sub) = gama_matrix(i_sub,j_sub,i_diag-1)
            end do 
        end do 
        aux = b_matrix - matmul(a_matrix,gama_aux)
        call inverse(aux,inv_matrix,dim)
        !
        !
        do i_sub = 1, dim
            d_v(i_sub) = res(i_sub,i_diag)
            beta_aux(i_sub) = beta(i_sub,i_diag-1)
        end do
        aux_v = d_v - matmul(a_matrix,beta_aux)
        !
        !
        beta_aux = matmul(inv_matrix,aux_v)
        !
        !
        do i_sub = 1, dim 
            beta(i_sub,i_diag) = beta_aux(i_sub)
        end do
        !
        !
    end do 
    !
    !
    i_diag = diag
    do i_sub = 1, dim 
        sol(i_sub,i_diag) = beta(i_sub,i_diag) 
    end do
    !
    !
    do i_diag = diag - 1, 1, -1
        !
        !
        do i_sub = 1, dim 
            x_v(i_sub) = sol(i_sub,i_diag+1)
        end do
        !
        !
        do j_sub = 1, dim
            do i_sub = 1, dim
                gama_aux(i_sub,j_sub) = gama_matrix(i_sub,j_sub,i_diag) 
            end do 
        end do 
        aux_v = matmul(gama_aux,x_v)
        !
        !
        do i_sub = 1, dim 
            sol(i_sub,i_diag) = beta(i_sub,i_diag)  - aux_v(i_sub)
        end do
        !
        !
    end do 
    !
    !
    ! do i_diag = 1, diag
    !     do i_sub = 1, dim
    !         write(*,*) sol(i_sub,i_diag)
    !     end do 
    ! end do 
    !
    !
end subroutine thomas_block_tridiagonal
!
!
!
subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer(4) :: n
double precision a(n,n), c(n,n), L(n,n), U(n,n)
double precision b(n),d(n),x(n)
!real(8),dimension(n,n) :: a, c, L, U
!real(8),dimension(n)   :: b,d,x
real(8)                :: coeff
integer(4)             :: i,j,k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do
! write by MarcoMatunaga
203 format (/, ' Lower Matrix')
204 format (/, ' Upper Matrix')
205 format (6f12.6)
! print the lower matrix L
  ! write(*,203) 
  ! do i=1,n
  !   write(*,205) (L(i,j),j=1,n)
  ! enddo
! print the upper matrix U
  ! write(*,204)
  ! do i=1,n
  !   write(*,205) (U(i,j),j=1,n)
  ! enddo
! end of MarcoMatunaga writing
! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
!
!
end subroutine inverse