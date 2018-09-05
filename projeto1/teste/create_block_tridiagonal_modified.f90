!
! retired from https://stackoverflow.com/questions/39256058/build-a-block-tri-diagonal-matrix
! write by High Performance Mark
! user information: https://stackoverflow.com/users/44309/high-performance-mark
!
program test
  implicit none 
  double precision lower(4,4), main(4,4), upper(4,4)
  integer, parameter :: jmax = 3
  integer, parameter :: n = 4
  integer i 

  data (lower(1,i), i=1,n) /  2.0,  1.0, 1.0, 0.0 /
  data (lower(2,i), i=1,n) /  4.0,  3.0, 3.0, 1.0 /
  data (lower(3,i), i=1,n) /  8.0,  7.0, 9.0, 5.0 /
  data (lower(4,i), i=1,n) /  6.0,  7.0, 9.0, 8.0 /

  data (main(1,i), i=1,n) /  2.2,  1.2, 1.2, 0.2 /
  data (main(2,i), i=1,n) /  4.2,  3.2, 3.2, 1.2 /
  data (main(3,i), i=1,n) /  8.2,  7.2, 9.2, 5.2 /
  data (main(4,i), i=1,n) /  6.2,  7.2, 9.2, 8.2 /

  data (upper(1,i), i=1,n) /  2.4,  1.4, 1.4, 0.4 /
  data (upper(2,i), i=1,n) /  4.4,  3.4, 3.4, 1.4 /
  data (upper(3,i), i=1,n) /  8.4,  7.4, 9.4, 5.4 /
  data (upper(4,i), i=1,n) /  6.4,  7.4, 9.4, 8.4 /

  call create_block_tridiagonal(lower,main,upper,4,jmax)
end program test
!
!
!
subroutine create_block_tridiagonal(l_diag,m_diag,u_diag,k,n)
    implicit none
  integer :: k  ! submatrix size
  integer :: n  ! number of submatrices along main diagonal
  INTEGER(4) :: ix, i, j   ! loop index
  integer :: id_column, id_line, id_block
  !
  ! we received the matrices as size x size
  double precision l_diag(k,k), m_diag(k,k), u_diag(k,k)
  ! for the block construction
  ! the submatrices, a (lower diagonal) b (main diagonal) c (upper diagonal)
  double precision amx(n-1,k,k), bmx(n,k,k), cmx(n-1,k,k)
  ! the block tridiagonal matrix
  double precision mat_a(n*k,n*k)
  !
  ! id_block = indixex of the block f
  ! for lower and upper diagonals is n - 1
  !
  mat_a = 0.0
  do id_block = 1, n - 1
      do id_line = 1, k
        do id_column = 1, k
            amx(id_block,id_line,id_column) = l_diag(id_line,id_column) 
            cmx(id_block,id_line,id_column) = u_diag(id_line,id_column) 
        end do 
      end do
  end do
  do id_block = 1, n
      do id_line = 1, k
        do id_column = 1, k
            bmx(id_block,id_line,id_column) = m_diag(id_line,id_column) 
        end do 
      end do
  end do
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


end subroutine create_block_tridiagonal
