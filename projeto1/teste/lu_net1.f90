program main
!====================================================================
!  Computing Inverse matrix
!  Method: Based on the Doolittle LU method
!====================================================================
implicit none
integer, parameter :: n=4
double precision a(n,n), c(n,n), U(n,n), L(n,n)
integer i,j
! matrix A
  data (a(1,i), i=1,n) /  2.0,  1.0, 1.0, 0.0 /
  data (a(2,i), i=1,n) /  4.0,  3.0, 3.0, 1.0 /
  data (a(3,i), i=1,n) /  8.0,  7.0, 9.0, 5.0 /
  data (a(4,i), i=1,n) /  6.0,  7.0, 9.0, 8.0 /

! print a header and the original matrix
  write (*,200)
  do i = 1 , n
     write (*,201) (a(i,j),j=1,n)
  end do
  
  call inverse(a,c,n)

! print the inverse matrix C = A^{-1} 
  write (*,202)
  do i = 1,n
     write (*,201)  (c(i,j),j=1,n)
  end do
200 format (' Computing Inverse matrix ',/,/, &
            ' Matrix A')
201 format (6f12.6)
202 format (/,' Inverse matrix A^{-1}')

end

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
real(8),dimension(n,n) :: a, c, L, U
real(8),dimension(n)   :: b,d,x
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
  write(*,203) 
  do i=1,n
    write(*,205) (L(i,j),j=1,n)
  enddo
! print the upper matrix U
  write(*,204)
  do i=1,n
    write(*,205) (U(i,j),j=1,n)
  enddo
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
end subroutine inverse