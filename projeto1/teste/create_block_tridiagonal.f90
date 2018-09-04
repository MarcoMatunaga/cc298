PROGRAM test
  IMPLICIT NONE

  INTEGER(4) :: k    ! submatrix size
  INTEGER(4) :: n    ! number of submatrices along main diagonal
  INTEGER(4) :: ix   ! loop index

  ! the submatrices, a (lower diagonal) b (main diagonal) c (upper diagonal)
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: amx, bmx, cmx

  ! the block tridiagonal matrix
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: mat_a

  k = 3  ! set these values as you wish
  n = 4

  ALLOCATE(amx(n-1,k,k), bmx(n,k,k), cmx(n-1,k,k))
  ALLOCATE(mat_a(n*k,n*k))

  mat_a = 0.0

  ! populate these as you wish
  amx = 1.0
  bmx = 2.0
  cmx = 3.0

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

open(1, file ='block_matrix')
write(1,*) mat_a
close(1)
END PROGRAM test

