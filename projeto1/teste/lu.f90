program lu
        implicit none
        real(8), dimension(4,4) :: a, alfa, beta
        real(8)                 :: sum, sum1
        integer(4)              :: i,j,k,n
!
! A - matrix to be decomposed
!
a(1,1) = 2.0d0
a(2,1) = 4.0d0
a(3,1) = 8.0d0
a(4,1) = 6.0d0
a(1,2) = 1.0d0
a(2,2) = 3.0d0
a(3,2) = 7.0d0
a(4,2) = 7.0d0
a(1,3) = 1.0d0
a(2,3) = 3.0d0
a(3,3) = 9.0d0
a(4,3) = 9.0d0
a(1,4) = 0.0d0
a(2,4) = 1.0d0
a(3,4) = 5.0d0
a(4,4) = 8.0d0
!
! Numerical recipes in FORTRAN77
! Chapter 2 - Solution of Linear Algebraic Equations
!
! algorithm of crout - decompose a matrix in two
! we use this as to create our upper and lower matrix
! L = lower matrix and U = upper matrix
!
! | beta(1,1) alfa(1,2) beta(1,3) beta(1,4) |
! | alfa(2,1) beta(2,2) beta(2,3) beta(2,4) |
! | alfa(3,1) alfa(3,2) beta(3,3) beta(3,4) |
! | alfa(4,1) alfa(4,2) alfa(4,3) beta(4,4) |
!
alfa = 0.0d0
beta = 0.0d0
sum  = 0.0d0
sum1 = 0.0d0 
n = 4 ! size of the matrix
!
! step 1 - 
!
do j = 1, n
    do i = 1, n
        if(i==j) alfa(i,j) = 1.0d0
    end do
end do
!
! step 2 -
!
do j = 1, n
    !
    ! 
        do i = 1, j 
            do k = 1, i - 1
                sum = sum + alfa(i,k)*beta(k,j)
            end do
            beta(i,j) = a(i,j) - sum 
        end do
    !
    !

        do i = j+1, n 
            do k = 1, j - 1
                sum1  = sum1 + alfa(i,k)*beta(k,j) 
            end do 
            alfa(i,j) = (1.0d0/beta(j,j))*(a(i,j) - sum1)
        end do
    !
    !
end do
!
!
print *, '|', beta(1,1), alfa(1,2), beta(1,3), beta(1,4), '|'
print *, '|', alfa(2,1), beta(2,2), beta(2,3), beta(2,4), '|'
print *, '|', alfa(3,1), alfa(3,2), beta(3,3), beta(3,4), '|'
print *, '|', alfa(4,1), alfa(4,2), alfa(4,3), beta(4,4), '|'
end program lu
