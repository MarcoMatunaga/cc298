!
! retired from https://stackoverflow.com/questions/39256058/build-a-block-tri-diagonal-matrix
! write by High Performance Mark
! modified by MarcoMatunaga
! user information: https://stackoverflow.com/users/44309/high-performance-mark
!
subroutine create_block_tridiagonal(l_diag,m_diag,u_diag,k,n)
    implicit none
  integer k  ! submatrix size
  integer n  ! number of submatrices along main diagonal
  integer ix, i, j   ! loop index
  integer id_column, id_line, id_block
  !
  ! we received the matrices as size x size 
  ! i.e., this subroutine transforms a quadratic matrix 
  ! into a block tridiagonal matrix
  !
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
  !
  !
  !
  do id_block = 1, n - 1
      do id_line = 1, k
        do id_column = 1, k
            amx(id_block,id_line,id_column) = l_diag(id_line,id_column) 
            cmx(id_block,id_line,id_column) = u_diag(id_line,id_column) 
        end do 
      end do
  end do
  !
  !
  !
  do id_block = 1, n
      do id_line = 1, k
        do id_column = 1, k
            bmx(id_block,id_line,id_column) = m_diag(id_line,id_column) 
        end do 
      end do
  end do
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
! 201 format (20f12.6)
!
!
end subroutine create_block_tridiagonal