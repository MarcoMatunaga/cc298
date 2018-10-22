subroutine pulliam_chausse
    use vars
    use diagonalization
    implicit none
    !integer(4)                                     :: i_sol, j_sol
    real(8)                                        :: rho_t
    real(8),dimension(:,:,:),allocatable           :: Id_x, Id_y
    real(8),dimension(:,:),allocatable             :: inverse_T_ksi
    real(8),dimension(:,:),allocatable             :: lambda_plus, lambda_minus
    real(8),dimension(:,:),allocatable             :: Identy
    real(8),dimension(:),allocatable               :: aux_residue
    real(8),dimension(:),allocatable               :: result
!
!
allocate(Identy(dim,dim))
allocate(inverse_T_ksi(dim,dim))
allocate(lambda_minus(dim,dim),lambda_plus(dim,dim))
allocate(Id_x(dim,dim,imax),Id_y(dim,dim,jmax))
allocate(aux_residue(dim))
allocate(result(dim))
!
!
do i = 1, dim
     Id_x(i,i,1:imax) = 1.0d0 
     Id_y(i,i,1:jmax) = 1.0d0
     Identy(i,i) = 1.0d0 
end do
!
!
j = 2
!
!
    i = 2   
    rho_t         = Q_barra(i+1,j,1)/metric_jacobian(i+1,j)
    inverse_T_ksi = inv_T_ksi(u(i,j),v(i,j),rho_t,a(i,j),ksi_x(i,j),ksi_y(i,j),dim)  
        call compute_residue(i,j)
        aux_residue(1) = residue(i,j,1)
        aux_residue(2) = residue(i,j,2)
        aux_residue(3) = residue(i,j,3)
        aux_residue(4) = residue(i,j,4)
    result = matmul(inverse_T_ksi,aux_residue)
    !
    !  
!
!
deallocate(Identy)
deallocate(inverse_T_ksi)
deallocate(lambda_minus,lambda_plus)
deallocate(Id_x,Id_y)
deallocate(aux_residue)
deallocate(result)
!
!
end subroutine pulliam_chausse