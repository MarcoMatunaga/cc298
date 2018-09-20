subroutine thomas_block_tridiagonal_v0(x_dir,y_dir,size,lower,main,upper,sol,res)
    use vars
    implicit none
!
! the system that will be solved is denoted by A X = B 
! where A is the block tridiagonal, x the vector of 'solutions', and B the vector
! of 'results'
!
! the lower diagonal is denote by a_matrix, 
! the main diagonal is denoted by b_matrix,
! the upper diagonal is denoted by c_matrix.
!
! lower, main and upper are matrices
!
! sol and res are matrices as well just to make
! easier the implementation (see pulliam to understand)
!
! size is the number of elemensts in the main diagonal
!
integer size, i_sub, j_sub, diag, num_call
integer x_dir, y_dir, index_dir
double precision sol(imax,jmax,dim), res(imax,jmax,dim)
double precision lower(size,dim,dim), main(size,dim,dim), upper(size,dim,dim)
!
! declare the variables used to solve the block tridiagonal system
! this algorithm follows the pdf in folder project1 (section4.pdf)
!
double precision b_matrix(dim,dim), inv_matrix(dim,dim), c_matrix(dim,dim), a_matrix(dim,dim)
double precision gama_matrix(size,dim,dim), gama_matrix_aux(dim,dim)
!
! vectors of the linear system
!
double precision beta(size,dim), x(size,dim), d(size,dim)
double precision beta_aux(dim), x_aux(dim), d_aux(dim)
!
! let is start
!
b_matrix = 0.0d0
c_matrix = 0.0d0
inv_matrix = 0.0d0
a_matrix = 0.0d0
gama_matrix = 0.0d0
beta_aux = 0.0d0
d_aux = 0.0d0
x_aux = 0.0d0
!
! gama(1) = (inverse of b1) * (c1)
!

    b_matrix(i_sub,j_sub) = main(1,i_sub,j_sub)

! gama_matrix(1,i_sub,j_sub) = 
!
! gama(k) = ( inverse of {b(k) - a(k)*gama(k-1)} ) * c(k)
! k = 2, n - 1
!

!
! beta(1) = ( inverse of b1) * (d1)
!

!
! beta(k) = ( inverse{b(k) - a(k)*gama(k-1)} )* ( d(k) - A(k)*beta(k-1) )
! k = 2, n 
!

!
! back substitution
! 
! x(n) = beta(n) 
!

!
! x(k) = beta(k) - gama(k)*x(k-1)
! k = (n-1), (n-2), 1
!
end subroutine thomas_block_tridiagonal_v0
