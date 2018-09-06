!
! retired from https://stackoverflow.com/questions/39256058/build-a-block-tri-diagonal-matrix
! write by High Performance Mark
! modified by MarcoMatunaga
! user information: https://stackoverflow.com/users/44309/high-performance-mark
!
subroutine create_block_tridiagonal(amx,bmx,cmx,k,n,mat_a)
    implicit none
  ! the block tridiagonal matrix
  double precision mat_a(n*k,n*k)
  integer k  ! submatrix size
  integer n  ! number of submatrices along main diagonal
  integer ix, i, j   ! loop index
  ! for the block construction
  ! the submatrices, a (lower diagonal) b (main diagonal) c (upper diagonal)
  double precision amx(n-1,k,k), bmx(n,k,k), cmx(n-1,k,k)
  !
  ! first the lower diagonal
  !
  DO ix = 1,k*(n-1),k
     mat_a(ix+k:ix+2*k-1,ix:ix+k-1) = amx(CEILING(REAL(ix)/REAL(k)),:,:)
  END DO
  !
  ! now the main diagonal
  !
  DO ix = 1,k*n,k
     mat_a(ix:ix+k-1,ix:ix+k-1) = bmx(CEILING(REAL(ix)/REAL(k)),:,:)
  END DO
  !
  ! finally the upper diagonal
  !
  DO ix = 1,k*(n-1),k
     mat_a(ix:ix+k-1,ix+k:ix+2*k-1) = cmx(CEILING(REAL(ix)/REAL(k)),:,:)
  END DO
!
! use for test or debug 
!
!   write (*,200)
!   do i = 1,n*k
!      write (*,201)  (mat_a(i,j),j=1,n*k)
!   end do
! 200 format (' Matrix A - block tridiagonal')
! !201 format (20f12.6)
!   write (*,200)
!   do i = 1,n*k
!      write (*,201) mat_a(i,i)
!   end do
! 200 format (' Matrix A - block tridiagonal')
! 201 format (20f12.6)
!
!
end subroutine create_block_tridiagonal