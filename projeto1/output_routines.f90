module output_routines
!
!
!
contains 
!
!
!        
        character (len=2) function convert(a)
        !
        implicit none
        !
        integer(4) :: a
        !
        if (a < 10) then
        convert = char(48+a)
        else if (a >= 10 .and. a < 20) then 
        convert = char(48+1)//char(48+mod(a,10))
        !
        !
        else if (a >= 20 .and. a < 30) then
        convert = char(48+2)//char(48+mod(a,10))
        !
        !
        else if (a >= 30 .and. a < 40) then
        convert = char(48+3)//char(48+mod(a,10))
        !
        !
        else if (a >= 40 .and. a < 50) then
        convert = char(48+4)//char(48+mod(a,10))
        !
        !
        else if (a >= 50 .and. a < 60) then
        convert = char(48+5)//char(48+mod(a,10))
        !
        !
        else if (a >= 60 .and. a < 70) then
        convert = char(48+6)//char(48+mod(a,10))
        !
        !
        else if (a >= 70 .and. a < 80) then
        convert = char(48+7)//char(48+mod(a,10))
        !
        !
        else if (a >= 80 .and. a < 90) then
        convert = char(48+8)//char(48+mod(a,10))
        !
        !
        else if (a >= 90 .and. a < 100) then
        convert = char(48+9)//char(48+mod(a,10))
        !
        !
        else
        write(*,*) 'Number os solutions saved more than 99'
        end if
        !
        !
        end function convert
    !
    !
    !
    subroutine output_metric_terms 
        use vars
        implicit none
    !
    !
    !
    open(6, file="metric_terms.dat")
    write(6,*) 'TITLE = "Projeto1" '
    write(6,*) 'VARIABLES = "x" "y" "i" "j" "x_ksi" "y_ksi" "x_eta" "y_eta" "metric_jacobian" ' 
    write(6,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT'
    do j = 1, jmax
        do i = 1, imax
            write(6,*) meshx(i,j), meshy(i,j), i,j, x_ksi(i,j), y_ksi(i,j), &
                                   x_eta(i,j), y_eta(i,j), metric_jacobian(i,j)
        end do 
    end do
    !
    !
    !
    close(6)
    end subroutine output_metric_terms
    !
    ! plot the fluxes E and F 
    !
    ! testar os fluxos
    !
    subroutine output_fluxes
        use vars
        implicit none
        !
        !
        !
    open(4, file='teste_fluxes.dat')
    write(4,*) 'TITLE = "Projeto1" '
    write(4,*) 'VARIABLES = "E_1" "E_2" "E_3" "E_4" '
    write(4,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
    do j = 1, jmax
        do i = 1, imax
            write(4,'(4EN20.10)') E_barra(i,j,1), E_barra(i,j,2), E_barra(i,j,3), E_barra(i,j,4) 
        end do
    end do
    close(4)
    !
    !
    !
    end subroutine output_fluxes
    !
    !
    !
    !
    subroutine output_inicial
        use vars 
        implicit none
        real(8),dimension(:,:),allocatable           :: p_out, u_out, v_out, rho, mach
    !
    !
    allocate(p_out(imax,jmax), u_out(imax,jmax), v_out(imax,jmax), rho(imax,jmax), mach(imax,jmax) )
    !
     do j = 1, jmax
         do i = 1, imax
             rho(i,j)   = Q_barra(i,j,1)/metric_jacobian(i,j)
             u_out(i,j) = Q_barra(i,j,2)/Q_barra(i,j,1)
             v_out(i,j) = Q_barra(i,j,3)/Q_barra(i,j,1)
             p_out(i,j) = (gama - 1.0d0)*( Q_barra(i,j,4)/metric_jacobian(i,j) &
                          - 0.5d0*rho(i,j)*( u_out(i,j)**2.0d0 + v_out(i,j)**2.0d0 ) )
             mach(i,j)  = sqrt( u_out(i,j)**2.0d0 + v_out(i,j)**2.0d0 ) / sqrt( gama*p_out(i,j)/rho(i,j) )
         end do
     end do
    !
    ! condicoes de contorno e iniciais
    !
    open(3,file='teste_init.dat')
    write(3,*) 'TITLE = "Projeto1" '
    write(3,*) 'VARIABLES =  "X" "Y" "i" "j" "p_curv" "p" "mach" '
    write(3,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
    do j = 1, jmax
        do i = 1, imax
            !write(3,'(7es11.3e2)') meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
            !write(3,*) meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
            write(3,'(9ES20.10)') meshx(i,j), meshy(i,j), DBLE(i), DBLE(j), p_out(i,j), p(i,j), mach(i,j)
        end do
    end do
    close(3)
    !
    !
    !
    deallocate(u_out, v_out, p_out, rho, mach)
    !
    !
    !
    end subroutine output_inicial
    !
    !
    !
    subroutine output_final
        use vars 
        implicit none
        real(8),dimension(:,:),allocatable           :: p_out, u_out, v_out, rho, mach
    !
    !
    allocate(p_out(imax,jmax), u_out(imax,jmax), v_out(imax,jmax), rho(imax,jmax), mach(imax,jmax) )
    !
     do j = 1, jmax
         do i = 1, imax
             rho(i,j)   = Q_barra(i,j,1)/metric_jacobian(i,j)
             u_out(i,j) = Q_barra(i,j,2)/Q_barra(i,j,1)
             v_out(i,j) = Q_barra(i,j,3)/Q_barra(i,j,1)
             p_out(i,j) = (gama - 1.0d0)*( Q_barra(i,j,4)/metric_jacobian(i,j) &
                          - 0.5d0*rho(i,j)*( u_out(i,j)**2.0d0 + v_out(i,j)**2.0d0 ) )
             mach(i,j)  = sqrt( u_out(i,j)**2.0d0 + v_out(i,j)**2.0d0 ) / sqrt( gama*p_out(i,j)/rho(i,j) )
         end do
     end do
    !
    ! condicoes de contorno e iniciais
    !
    open(3,file='teste_boundary.dat')
    write(3,*) 'TITLE = "Projeto1" '
    write(3,*) 'VARIABLES = "X" "Y" "u" "v" "rho" "p" "e" "mach" '
    write(3,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
    do j = 1, jmax
        do i = 1, imax
            !write(3,'(7es11.3e2)') meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
            !write(3,*) meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
            write(3,'(10ES20.10)') meshx(i,j), meshy(i,j), u_out(i,j), v_out(i,j), rho(i,j), p_out(i,j), & 
            Q_barra(i,j,4)/metric_jacobian(i,j), mach(i,j)
        end do
    end do
    close(3)
    !
    !
    !
    deallocate(u_out, v_out, p_out, rho, mach)
    !
    !
    !
    end subroutine output_final

    subroutine output_residue
        use vars
        implicit none
        !
        !
        !
    open(5, file="residue.dat")
    write(5,*) iter, max_residue
    end subroutine output_residue
!
!
!
    subroutine output_tecplot
        use vars 
        implicit none
        real(8),dimension(:,:),allocatable           :: p_out, u_out, v_out, rho, mach, q_vel_out
        character (len=100) :: fname = '.dat'
        character (len=100) :: FileTag    
    !
    !
    FileTag = convert(nsave)
    allocate(p_out(imax,jmax), u_out(imax,jmax), v_out(imax,jmax), rho(imax,jmax), mach(imax,jmax) )
    allocate(q_vel_out(imax,jmax))
    !
     do j = 1, jmax
         do i = 1, imax
             rho(i,j)   = Q_barra(i,j,1)/metric_jacobian(i,j)
             u_out(i,j) = Q_barra(i,j,2)/Q_barra(i,j,1)
             v_out(i,j) = Q_barra(i,j,3)/Q_barra(i,j,1)
             p_out(i,j) = (gama - 1.0d0)*( Q_barra(i,j,4)/metric_jacobian(i,j) &
                          - 0.5d0*rho(i,j)*( u_out(i,j)**2.0d0 + v_out(i,j)**2.0d0 ) )
             mach(i,j)  = sqrt( u_out(i,j)**2.0d0 + v_out(i,j)**2.0d0 ) / sqrt( gama*p_out(i,j)/rho(i,j) )
             q_vel_out(i,j) = sqrt( u_out(i,j)**2.0d0 + v_out(i,j)**2.0d0 )
         end do
     end do
    !
    ! condicoes de contorno e iniciais
    !
    open(7,file=trim(FileTag)//trim(fname))
    write(7,*) 'TITLE = "Projeto1" '
    write(7,*) 'VARIABLES = "X" "Y" "u" "v" "q_vel_out" "rho" "p" "e" "mach" '
    write(7,*) 'ZONE I = ', imax, ' J =', jmax, ' DATAPACKING = POINT' 
    do j = 1, jmax
        do i = 1, imax
            !write(3,'(7es11.3e2)') meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
            !write(3,*) meshx(i,j), meshy(i,j), x_ksi(i,j), x_eta(i,j), y_ksi(i,j), y_eta(i,j), metric_jacobian(i,j)
            write(7,'(15ES20.10)') meshx(i,j), meshy(i,j), u_out(i,j), v_out(i,j), q_vel_out(i,j), rho(i,j), p_out(i,j), & 
            Q_barra(i,j,4)/metric_jacobian(i,j), mach(i,j)
        end do
    end do
    close(7)
    !
    !
    !
    deallocate(u_out, v_out, p_out, rho, mach)
    deallocate(q_vel_out)
    !
    !
    end subroutine
    !
    !
    !
end module output_routines