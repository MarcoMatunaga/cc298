subroutine thomas_pulliam_chausse(a,b,c,d,x,n)
    implicit none
    !    a - sub-diagonal (means it is the diagonal below the main diagonal)
    !    b - the main diagonal
    !    c - sup-diagonal (means it is the diagonal above the main diagonal)
    !    d - right part
    !    x - the answer
    !    n - number of equations

        integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: a,b,c,d
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: cp,dp
        real(8) :: m
        integer tc_i
    ! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
    ! solve for vectors c-prime and d-prime
         do tc_i = 2,n
           m = b(tc_i)-cp(tc_i-1)*a(tc_i)
           cp(tc_i) = c(tc_i)/m
           dp(tc_i) = (d(tc_i)-dp(tc_i-1)*a(tc_i))/m
         enddo
    ! initialize x
         x(n) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
        do tc_i = n-1, 1, -1
          x(tc_i) = dp(tc_i)-cp(tc_i)*x(tc_i+1)
        end do
end subroutine thomas_pulliam_chausse