subroutine inv(A,A_inv,m)
  implicit none
  integer :: m
  real(8), dimension(m,m)::A, A_inv
  real(8),dimension(m)::WORK
  integer,dimension(m)::IPIV
  integer info
  A_inv = A
  call DGETRF(M,M,A_inv,M,IPIV,info)
  if (info /=  0) then
    write(*,*)"DGETRF: Failed during matrix factorization"
    stop
  end if
  call DGETRI(M,A_inv,M,IPIV,WORK,M,info)
  if (info /=  0) then
   write(*,*)"DGETRI: Failed during matrix inversion."
   stop
  end if
end subroutine inv