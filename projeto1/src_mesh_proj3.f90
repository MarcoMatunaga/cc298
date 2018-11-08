subroutine mesh_proj3
    use vars
    use vars_proj3
    implicit none
    real(8)                      :: Deltax, Deltay
    real(8),dimension(imax)      :: aux_x
    real(8),dimension(jmax)      :: aux_y
    integer(4)                   :: msh_i, msh_j

    Deltax = Length/DBLE(imax)
    Deltay = Height/DBLE(jmax)

    do msh_i = 1, imax
        aux_x(msh_i) = msh_i*Deltax
    end do

    do msh_j = 1, jmax
        aux_y(msh_j) = msh_j*Deltay
    end do

    do msh_j = 1, jmax
        do msh_i = 1, imax
            meshx(msh_i,msh_j) = aux_x(msh_i)
            meshy(msh_i,msh_j) = aux_y(msh_j)
        end do
    end do

    ! arquivo tecplot

    open(2,file='mesh.dat')
    write(2,*) 'TITLE = "Projeto1" '
    write(2,*) 'VARIABLES = "X" "Y" '
    write(2,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
    do msh_j = 1, jmax
        do msh_i = 1, imax
            write(2,*) meshx(msh_i,msh_j), meshy(msh_i,msh_j)
        end do
    end do

end subroutine mesh_proj3