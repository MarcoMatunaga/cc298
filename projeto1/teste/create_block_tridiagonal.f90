PROGRAM test
  IMPLICIT NONE

  integer, parameter :: k = 4 ! submatrix size
  integer, parameter :: n = 4 ! number of submatrices along main diagonal
  INTEGER(4) :: ix, i, j   ! loop index

  ! the submatrices, a (lower diagonal) b (main diagonal) c (upper diagonal)
!  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: amx, bmx, cmx
  double precision amx(n-1,k,k), bmx(n,k,k), cmx(n-1,k,k)
  ! the block tridiagonal matrix

!  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: mat_a
  double precision mat_a(n*k,n*k)
  ! k = 4  ! set these values as you wish
  ! n = 4

  ! ALLOCATE(amx(n-1,k,k), bmx(n,k,k), cmx(n-1,k,k))
  ! ALLOCATE(mat_a(n*k,n*k))

  mat_a = 0.0

  ! populate these as you wish

  data (amx(1,1,i), i=1,k) /  2.0,  1.0, 1.0, 0.0 /
  data (amx(1,2,i), i=1,k) /  4.0,  3.0, 3.0, 1.0 /
  data (amx(1,3,i), i=1,k) /  8.0,  7.0, 9.0, 5.0 /
  data (amx(1,4,i), i=1,k) /  6.0,  7.0, 9.0, 8.0 /

  data (cmx(1,1,i), i=1,k) /  2.4,  1.4, 1.4, 0.4 /
  data (cmx(1,2,i), i=1,k) /  4.4,  3.4, 3.4, 1.4 /
  data (cmx(1,3,i), i=1,k) /  8.4,  7.4, 9.4, 5.4 /
  data (cmx(1,4,i), i=1,k) /  6.4,  7.4, 9.4, 8.4 /

  data (bmx(1,1,i), i=1,k) /  2.2,  1.2, 1.2, 0.2 /
  data (bmx(1,2,i), i=1,k) /  4.2,  3.2, 3.2, 1.2 /
  data (bmx(1,3,i), i=1,k) /  8.2,  7.2, 9.2, 5.2 /	
  data (bmx(1,4,i), i=1,k) /  6.2,  7.2, 9.2, 8.2 /
!***************************************************!
  data (amx(2,1,i), i=1,k) /  2.0,  1.0, 1.0, 0.0 /
  data (amx(2,2,i), i=1,k) /  4.0,  3.0, 3.0, 1.0 /
  data (amx(2,3,i), i=1,k) /  8.0,  7.0, 9.0, 5.0 /
  data (amx(2,4,i), i=1,k) /  6.0,  7.0, 9.0, 8.0 /

  data (cmx(2,1,i), i=1,k) /  2.4,  1.4, 1.4, 0.4 /
  data (cmx(2,2,i), i=1,k) /  4.4,  3.4, 3.4, 1.4 /
  data (cmx(2,3,i), i=1,k) /  8.4,  7.4, 9.4, 5.4 /
  data (cmx(2,4,i), i=1,k) /  6.4,  7.4, 9.4, 8.4 /

  data (bmx(2,1,i), i=1,k) /  2.2,  1.2, 1.2, 0.2 /
  data (bmx(2,2,i), i=1,k) /  4.2,  3.2, 3.2, 1.2 /
  data (bmx(2,3,i), i=1,k) /  8.2,  7.2, 9.2, 5.2 /
  data (bmx(2,4,i), i=1,k) /  6.2,  7.2, 9.2, 8.2 /
  !***************************************************!
  data (amx(3,1,i), i=1,k) /  2.0,  1.0, 1.0, 0.0 /
  data (amx(3,2,i), i=1,k) /  4.0,  3.0, 3.0, 1.0 /
  data (amx(3,3,i), i=1,k) /  8.0,  7.0, 9.0, 5.0 /
  data (amx(3,4,i), i=1,k) /  6.0,  7.0, 9.0, 8.0 /

  data (cmx(3,1,i), i=1,k) /  2.4,  1.4, 1.4, 0.4 /
  data (cmx(3,2,i), i=1,k) /  4.4,  3.4, 3.4, 1.4 /
  data (cmx(3,3,i), i=1,k) /  8.4,  7.4, 9.4, 5.4 /
  data (cmx(3,4,i), i=1,k) /  6.4,  7.4, 9.4, 8.4 /

  data (bmx(3,1,i), i=1,k) /  2.2,  1.2, 1.2, 0.2 /
  data (bmx(3,2,i), i=1,k) /  4.2,  3.2, 3.2, 1.2 /
  data (bmx(3,3,i), i=1,k) /  8.2,  7.2, 9.2, 5.2 /
  data (bmx(3,4,i), i=1,k) /  6.2,  7.2, 9.2, 8.2 /
!***************************************************!
  data (bmx(4,1,i), i=1,k) /  2.2,  1.2, 1.2, 0.2 /
  data (bmx(4,2,i), i=1,k) /  4.2,  3.2, 3.2, 1.2 /
  data (bmx(4,3,i), i=1,k) /  8.2,  7.2, 9.2, 5.2 /
  data (bmx(4,4,i), i=1,k) /  6.2,  7.2, 9.2, 8.2 /
!***************************************************!

  ! first the lower diagonal
  DO ix = 1,k*(n-1),k
     mat_a(ix+k:ix+2*k-1,ix:ix+k-1) = amx(CEILING(REAL(ix)/REAL(k)),:,:)
  END DO

  ! now the main diagonal
  DO ix = 1,k*n,k
     mat_a(ix:ix+k-1,ix:ix+k-1) = bmx(CEILING(REAL(ix)/REAL(k)),:,:)
  END DO

  ! finally the upper diagonal
  DO ix = 1,k*(n-1),k
     mat_a(ix:ix+k-1,ix+k:ix+2*k-1) = cmx(CEILING(REAL(ix)/REAL(k)),:,:)
  END DO

  write (*,200)
  do i = 1,n*k
     write (*,201)  (mat_a(i,j),j=1,n*k)
  end do
200 format (' Matrix A - block tridiagonal')
201 format (20f12.6)

END PROGRAM test

