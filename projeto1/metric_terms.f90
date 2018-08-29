!
! metric terms
!
subroutine metric_terms
use vars
implicit none
!****************************************
! delta_eta e delta_ksi igual a 1: podemos colocar
! como delta_eta e delta_ksi como inputs
!****************************************
!
! calculo das derivadas nos pontos interiores da malha
!
do j = 2, jmax - 1
    do i = 2, imax - 1 
        x_eta(i,j) = ( meshx(i,j+1) - meshx(i,j-1) )/( 2.0d0*delta_eta )
        y_eta(i,j) = ( meshy(i,j+1) - meshy(i,j-1) )/( 2.0d0*delta_eta )
        x_ksi(i,j) = ( meshx(i+1,j) - meshx(i-1,j) )/( 2.0d0*delta_ksi )
        y_ksi(i,j) = ( meshy(i+1,j) - meshy(i-1,j) )/( 2.0d0*delta_ksi )
    end do
end do
!
! calculo das derivadas nas fronteiras usar diferancas centradas 
! para diminuir o erro de discretizacao espacial
! 
! lembrando que as fronteiras sao i,j = 1 e i,j = imax,jmax  
!
do j = 2, jmax - 1
    !
    ! calculo derivadas na direcao eta, onde y varia - diferencas centradas 
    !
    x_eta(imax,j) = ( meshx(imax,j+1) - meshx(imax,j-1) )/( 2.0d0*delta_eta )
    y_eta(imax,j) = ( meshy(imax,j+1) - meshy(imax,j-1) )/( 2.0d0*delta_eta )
    x_eta(1,j) = ( meshx(1,j+1) - meshx(1,j-1) )/( 2.0d0*delta_eta )
    y_eta(1,j) = ( meshy(1,j+1) - meshy(1,j-1) )/( 2.0d0*delta_eta )
    !
    ! calculo derivadas na direcao ksi, onde x varia - one sided
    !
    x_ksi(imax,j) = ( meshx(imax,j) - meshx(imax-1,j) )/( delta_ksi )
    y_ksi(imax,j) = ( meshy(imax,j) - meshy(imax-1,j) )/( delta_ksi )
    x_ksi(1,j) = ( meshx(2,j) - meshx(1,j) )/( delta_ksi )
    y_ksi(1,j) = ( meshy(2,j) - meshy(1,j) )/( delta_ksi )
end do
!
!
!
do i = 2, imax - 1
    !
    ! calculo das derivadas na direcao eta, onde y varia - one sided
    !
    x_eta(i,jmax) = ( meshx(i,jmax) - meshx(i,jmax-1) ) /( delta_eta )
    y_eta(i,jmax) = ( meshy(i,jmax) - meshy(i,jmax-1) ) /( delta_eta )
    x_eta(i,1) = ( meshx(i,2) - meshx(i,1) )/( delta_ksi )
    y_eta(i,1) = ( meshy(i,2) - meshy(i,1) )/( delta_ksi )
    !
    ! calculo das derivadas na direcao ksi, onde x varia - diferencas centradas
    !
    x_ksi(i,jmax) = ( meshx(i+1,jmax) - meshx(i-1,jmax) )/( 2.0d0*delta_eta )
    y_ksi(i,jmax) = ( meshy(i+1,jmax) - meshy(i-1,jmax) )/( 2.0d0*delta_eta )
    x_ksi(i,1) = ( meshx(i+1,1) - meshx(i-1,1) )/( 2.0d0*delta_eta )
    y_ksi(i,1) = ( meshy(i+1,1) - meshy(i-1,1) )/( 2.0d0*delta_eta )
end do
!
! ainda falta os quatro pontos nos cantos do dominio computacional
! (1,jmax) (imax,jmax) (1,1) (imax,1)
!
! calculo no ponto (1,jmax)
!
  x_eta(1,jmax) = ( meshx(1,jmax) - meshx(1,jmax-1) )/( delta_eta )
  x_ksi(1,jmax) = ( meshx(2,jmax) - meshx(1,jmax) )/( delta_ksi )
  y_eta(1,jmax) = ( meshy(1,jmax) - meshy(1,jmax-1) )/( delta_eta )
  y_ksi(1,jmax) = ( meshy(2,jmax) - meshy(1,jmax) )/( delta_ksi )
! calculo no ponto (imax,jmax)
!
  x_eta(imax,jmax) = ( meshx(imax,jmax) - meshx(imax,jmax-1) ) / ( delta_eta )
  x_ksi(imax,jmax) = ( meshx(imax,jmax) - meshx(imax-1,jmax) ) / ( delta_ksi )
  y_eta(imax,jmax) = ( meshy(imax,jmax) - meshy(imax,jmax-1) ) / ( delta_eta ) 
  y_ksi(imax,jmax) = ( meshy(imax,jmax) - meshy(imax-1,jmax) ) / ( delta_ksi )
! calculo no ponto (1,1)
!
  x_eta(1,1) = ( meshx(1,2) - meshx(1,1) ) / ( delta_eta )
  x_ksi(1,1) = ( meshx(2,1) - meshx(1,1) ) / ( delta_ksi )
  y_eta(1,1) = ( meshy(1,2) - meshy(1,1) ) / ( delta_eta )
  y_ksi(1,1) = ( meshy(2,1) - meshy(1,1) ) / ( delta_ksi )
! calculo no ponto (imax,1)
!
  x_eta(imax,1) = ( meshx(imax,2) - meshx(imax,1) ) / ( delta_eta )
  x_ksi(imax,1) = ( meshx(imax,1) - meshx(imax-1,1) ) / ( delta_ksi )
  y_eta(imax,1) = ( meshy(imax,2) - meshy(imax,1) ) / ( delta_eta )
  y_ksi(imax,1) = ( meshy(imax,1) - meshy(imax-1,1) ) / ( delta_ksi )
!
!
!
do j = 1, jmax
    do i = 1, imax
        !
        ! be carefyl the fluxes vectors are multiplied by the
        ! x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j) not by 
        ! 1/x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j)
        !
        ! metric_jacobian(i,j) = 1.0d0/(x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j))
        metric_jacobian(i,j) = x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j)
        eta_x(i,j) = -1.0d0*(1.0d0/metric_jacobian(i,j))*y_ksi(i,j)
        eta_y(i,j) = (1.0d0/metric_jacobian(i,j))*x_ksi(i,j)
        ksi_x(i,j) = (1.0d0/metric_jacobian(i,j))*y_eta(i,j)
        ksi_y(i,j) = -1.0d0*(1.0d0/metric_jacobian(i,j))*x_eta(i,j)
    end do
end do
!
!
end subroutine metric_terms