module functionsDerivatives
    use vars
    implicit none
    
    contains
    
     subroutine FDFluxJacobian(qbarra1, qbarra2, qbarra3, qbarra4, ksix, ksiy, step, Jac)
          implicit none
        
          real(8), intent(in) :: qbarra1, qbarra2, qbarra3, qbarra4, ksix, ksiy, step
          real(8), intent(out) ::  Jac(dim,dim)
        
          real(8)             :: f1, f2
          real(8), dimension(dim)             :: qbarra

          Jac = 0.0d0
          
          ! linha 1
          qbarra(1) = qbarra1 + step
          f1 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          qbarra(1) = qbarra1 - step
          f2 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          Jac(1,1) = (f1 - f2)/(2.0d0*step)
        
          qbarra(2) = qbarra2 + step        
          f1 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          qbarra(2) = qbarra2 - step
          f2 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          Jac(1,2) = (f1 - f2)/(2.0d0*step)
        
          qbarra(3) = qbarra3 + step
          f1 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          qbarra(3) = qbarra3 - step
          f2 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          Jac(1,3) = (f1 - f2)/(2.0d0*step)
        
          qbarra(4) = qbarra4 + step
          f1 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          qbarra(4) = qbarra4 - step
          f2 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          Jac(1,4) = (f1 - f2)/(2.0d0*step)

          ! linha 2
          qbarra(1) = qbarra1 + step
          f1 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          qbarra(1) = qbarra1 - step
          f2 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) & 
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          Jac(2,1) = (f1 - f2)/(2.0d0*step)
        
          qbarra(2) = qbarra2 + step        
          f1 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          qbarra(2) = qbarra2 - step
          f2 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) & 
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          Jac(2,2) = (f1 - f2)/(2.0d0*step)
        
          qbarra(3) = qbarra3 + step
          f1 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          qbarra(3) = qbarra3 - step
          f2 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) & 
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          Jac(2,3) = (f1 - f2)/(2.0d0*step)
        
          qbarra(4) = qbarra4 + step
          f1 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          qbarra(4) = qbarra4 - step
          f2 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) & 
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          Jac(2,4) = (f1 - f2)/(2.0d0*step)

          ! linha 3
          qbarra(1) = qbarra1 + step
          f1 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          qbarra(1) = qbarra1 - step
          f2 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) & 
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          Jac(3,1) = (f1 - f2)/(2.0d0*step)
        
          qbarra(2) = qbarra2 + step        
          f1 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          qbarra(2) = qbarra2 - step
          f2 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) & 
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          Jac(3,2) = (f1 - f2)/(2.0d0*step)
        
          qbarra(3) = qbarra3 + step
          f1 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          qbarra(3) = qbarra3 - step
          f2 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) & 
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          Jac(3,3) = (f1 - f2)/(2.0d0*step)
        
          qbarra(4) = qbarra4 + step
          f1 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          qbarra(4) = qbarra4 - step
          f2 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) & 
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          Jac(3,4) = (f1 - f2)/(2.0d0*step)

          ! linha 4
          qbarra(1) = qbarra1 + step
          f1 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )
          qbarra(1) = qbarra1 - step
          f2 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )
          Jac(4,1) = (f1 - f2)/(2.0d0*step)
        
          qbarra(2) = qbarra2 + step        
          f1 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )
          qbarra(2) = qbarra2 - step
          f2 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )
          Jac(4,2) = (f1 - f2)/(2.0d0*step)
        
          qbarra(3) = qbarra3 + step
          f1 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )             
          qbarra(3) = qbarra3 - step
          f2 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )             
          Jac(4,3) = (f1 - f2)/(2.0d0*step)
        
          qbarra(4) = qbarra4 + step
          f1 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )
          qbarra(4) = qbarra4 - step
          f2 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )
          Jac(4,4) = (f1 - f2)/(2.0d0*step)

     end subroutine FDFluxJacobian
    
     ! return np.imag(F(x + 1j*h))/h
     subroutine complexStepMethod(qbarra1, qbarra2, qbarra3, qbarra4, ksix, ksiy, step, Jac)
          implicit none
          real(8), intent(in) :: qbarra1, qbarra2, qbarra3, qbarra4, ksix, ksiy, step
          real(8), intent(out) ::  Jac(dim,dim)
          !
          complex(16)                             :: step_cplx, f1
          complex(16), dimension(dim)             :: qbarra

          step_cplx = complex(0.0d0, step)
          Jac = 0.0d0
          !
          qbarra(1) = qbarra1 + step_cplx
          f1 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          Jac(1,1) = (aimag(f1)/(step))
        
          qbarra(2) = qbarra2 + step_cplx        
          f1 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          Jac(1,2) = (aimag(f1)/(step))
        
          qbarra(3) = qbarra3 + step_cplx
          f1 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          Jac(1,3) = (aimag(f1)/(step))
        
          qbarra(4) = qbarra4 + step_cplx
          f1 = (qbarra(2)*ksix + qbarra(3)*ksiy)
          Jac(1,4) = (aimag(f1)/(step))

          ! linha 2
          qbarra(1) = qbarra1 + step_cplx
          f1 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          Jac(2,1) = (aimag(f1)/(step))
        
          qbarra(2) = qbarra2 + step_cplx        
          f1 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          Jac(2,2) = (aimag(f1)/(step))
        
          qbarra(3) = qbarra3 + step_cplx
          f1 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          Jac(2,3) = (aimag(f1)/(step))
        
          qbarra(4) = qbarra4 + step_cplx
          f1 = (qbarra(2)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksix*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          Jac(2,4) = (aimag(f1)/(step))

          ! linha 3
          qbarra(1) = qbarra1 + step_cplx
          f1 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ) )
          Jac(3,1) = (aimag(f1)/(step))
        
          qbarra(2) = qbarra2 + step_cplx        
          f1 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          Jac(3,2) = (aimag(f1)/(step))
        
          qbarra(3) = qbarra3 + step_cplx
          f1 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          Jac(3,3) = (aimag(f1)/(step))
        
          qbarra(4) = qbarra4 + step_cplx
          f1 = (qbarra(3)/qbarra(1))*(qbarra(2)*ksix + qbarra(3)*ksiy) &
               + ksiy*((gama - 1.0d0)*(qbarra(4) & 
               - 0.50d0*qbarra(1)*(qbarra(2)**2 + qbarra(3)**2)/qbarra(1)**2 ))
          Jac(3,4) = (aimag(f1)/(step))

          ! linha 4
          qbarra(1) = qbarra1 + step_cplx
          f1 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )
          Jac(4,1) = (aimag(f1)/(step))
        
          qbarra(2) = qbarra2 + step_cplx        
          f1 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )
          Jac(4,2) = (aimag(f1)/(step))
        
          qbarra(3) = qbarra3 + step_cplx
          f1 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )             
          Jac(4,3) = (aimag(f1)/(step))
        
          qbarra(4) = qbarra4 + step_cplx
          f1 = ( (qbarra(2)/qbarra(1))*ksix + (qbarra(3)/qbarra(1))*ksiy )* &
               (qbarra(4) + (gama-1.0d0)*(qbarra(4) - 0.50d0*(1.0d0/qbarra(1))* &
               qbarra(2)**2 + qbarra(3)**2 ) )
          Jac(4,4) = (aimag(f1)/(step))

    end subroutine complexStepMethod

end module functionsDerivatives