
SUBROUTINE Derivs(t, depv, ddt)
!--------------------------------------------------------
! first order ODE system to integrate
!    drsx/dt=vsx
!    dvsx/dt = (G*mp/rsep^3)*(rpx-rsx)
!    similar for y components and for second  planet
!    possibly add 1/r^3 perturbation
! Derivs is called by module equdif.f90
! ddt = d()/dt
! depv = (rsx,vsx,rsy,vsy,rpx,vpx,rpy,vpy)
!--------------------------------------------------------HSJ


USE com,                    ONLY : rsep,cgs,cgp
IMPLICIT NONE

REAL*8  depv(*),ddt(*),t
  rsep   = SQRT((depv(5)-depv(1))**2+ (depv(7)-depv(3))**2)
  ddt(1) = depv(2)  
  ddt(2) = (cgs/rsep**3)*(depv(5)-depv(1)) 
  ddt(3) = depv(4)
  ddt(4) = (cgs/rsep**3)*(depv(7)-depv(3))
  ddt(5) = depv(6)
  ddt(6) = (cgp/rsep**3)*(depv(1)-depv(5))
  ddt(7) = depv(8)
  ddt(8) = (cgp/rsep**3)*(depv(3)-depv(7))
!  PRINT *,'depv(2),ddt(2)=',depv(2),ddt(2)
RETURN
END SUBROUTINE Derivs
