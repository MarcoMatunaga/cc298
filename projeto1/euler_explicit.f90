!
! euler explicit
!
subroutine euler_explicit
        use vars
        implicit none
!
! ************************
! Q_dis eh usado apenas aqui, podemos mudar 
! a declaracao da variavel para esse local
! ************************
!
do j = 1, jmax
        do i = 1, imax
            Q_dis(i,j,1) = Q_barra(i,j,1)/metric_jacobian(i,j)
            Q_dis(i,j,2) = Q_barra(i,j,2)/metric_jacobian(i,j)
            Q_dis(i,j,3) = Q_barra(i,j,3)/metric_jacobian(i,j)
            Q_dis(i,j,4) = Q_barra(i,j,4)/metric_jacobian(i,j)
        end do
end do
!
!
max_residue = -100.0d0
!
! explicit time marching 
!
do j = 2, jmax - 1
        do i = 2, imax - 1
            call residue(i,j)
            Q_barra(i,j,1) = Q_barra(i,j,1) - residue1(i,j)
            Q_barra(i,j,2) = Q_barra(i,j,2) - residue2(i,j)
            Q_barra(i,j,3) = Q_barra(i,j,3) - residue3(i,j) 
            Q_barra(i,j,4) = Q_barra(i,j,4) - residue4(i,j) 
        end do
end do
!
!
!
end subroutine euler_explicit