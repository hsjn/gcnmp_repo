
PROGRAM Grav_sling
!-------------------------------------------------------------------------
!*    CALLS   SUBROUTINE ODEINT OF MODULE EQUDIF    
!*       (Runge-Kutta Method with Time Step Control)        
!* -----------------------------------------------------------------------     
!*                                                                        *
!* NOTE that the following derivatives                                    *
!*      are given in Subroutine Derivs(t, XT, XF)                         *
!*      which is stored in file desys.f90                                 *
!* Integrate First Order ODE System:                                      *
!*    dependent variables are                                             *
!*   (rsx,vsx,rsy,vsy,rpx,vpx,rpy,vpy)                                    *
!*    where s = satellite,p = planet                                      *
!*                                                                        *
!*   Y1' = Y(2)                         ! drsx/dt = vsx                   *
!*   Y2' = (cgs/rsep^3)*(Y(5)-y(1))     ! dvsx/dt = G*mp/r^3 dot unitx    *
!*   Y3' = Y(4)                         ! drsy/dt = vsy                   *
!*   Y4' = (cgs/rsep^3)*(Y(7)-Y(3)      ! dvsy/dt = G*mp/r^3 dot unity    *
!*   Y5' = Y(6)                         ! drpx/dt = vpx                   *
!*   Y6' = (cgp/rsep^3)(Y(1)-Y(5))      ! dvpx/dt = G*ms/r^3 dot -unitx   *
!*   Y7' = Y(8)                         ! drpy/dt = vpy                   *
!*   Y8' = (cgp/rsep^3)(Y(3)-Y(7))      ! dvpy/dt = G*ms/r^3 dot -unity   *
!*                                                                        *
!*   from t1=0 to t2=10 with additional data:                             *
!*   Y start = (rs0x,vs0x,rs0y,vs0y,rp0x,vp0x,rp0y,vp0y)                  *
!-------------------------------------------------------------------------HSJ

USE Equdif                                !  For subroutine odein
USE com,                       ONLY : cgs,cgp,rsep
IMPLICIT NONE
INTEGER, PARAMETER :: NEQ = 8  
INTEGER  nok,nbad                         ! # eq, # ok time steps,# bad time steps
INTEGER  i,j,nsteps
REAL*8   Y(NEQ)
REAL*8   dt,dtf,dts,t1,t2,tol,secpday,tend_sec,tend_day
REAL*8,  PARAMETER ::  G = 6.67408e-11             ! M^3 kg^-1 s^-2 =N m^2/kg^2
REAL*8,  PARAMETER ::  ms = 10000.                 ! mass sattelite,kg
REAL*8,  PARAMETER ::  me = 5.97e24                ! mass earth,kg
REAL*8   rs,rp,rsepmin,rsepmax,rsepx,rsepy,rsep0   ! r sat, r plan, r separation,m
REAL*8   rs0x,vs0x,rp0x,vp0x                       ! x,y planar motion
REAL*8   rsrelx,rsrely
REAL*8   vs,vp,vrel                                ! ditto  speed,m/s
REAL*8   rs0y,vs0y,rp0y,vp0y
REAL*8   engs,engp,engs0,engp0,engt                ! kinetic eng sat and planet
REAL*8   mp,angmtm                                 ! mass planet,ang mtm of system about origin
! initial conditions:
   rs0x      = .766044e9  ; rs0y = -.64279e9
   vs0x      = -.766044e4 ; vs0y =  .64279e4
   vs0x      = -.766044e3 ; vs0y =  .64279e4


   engs0     = 0.5D0*ms*(vs0x**2+vs0y**2)

   mp        = 100.*me      ! me is earth mass kg
   mp        = 1.9e27       ! mass Jupiter
   rp0x      = 0.0  ; rp0y =  0.0
   vp0x      = 1.3e4 ; vp0y =  0.0 
   engp0     = 0.5D0*mp*(vp0x**2+vp0y**2)


!  rsep points from sat to planet
   rsep      = SQRT((rp0x-rs0x)**2 + (rp0y - rs0y)**2)
   rsepmin   = rsep ; rsepmax =0.0 ; rsep0 = rsep
   rsepx     = rp0x -rs0x
   rsepy     = rp0y -rs0y    
   cgp       = G*ms ; cgs = G*mp
   engt      = engp0 + engs0
   angmtm    = ms*(rs0x*vs0y - rs0y*vs0x) + mp*(rp0x*vp0y - rp0y*vp0x) 
   
   nsteps    = 0
   nok       =0   ; nbad   = 0



   OPEN(unit=9,file="out.txt",status="unknown")
!
!
  Y(1) = rs0x ; Y(2) = vs0x ; Y(3) = rs0y ; y(4)= vs0y
  Y(5) = rp0x ; Y(6) = vp0x ; Y(7) = rp0y ; y(8)= vp0y


  secpday  = 86400.   ! sec p day
  tend_day = 100.
  tend_sec = tend_day*secpday
  dt=30.d0            !starting integration step ,sec
  dts = dt
  t1=0.d0              !starting time
  t2=t1+1000.*dt       !ending time  
  tol = 1.0d-8


!----------------------------------------------------------------
!  Subroutine odeint(ystart, nvar, t1, t2, eps, h1, hmin, nok, nbad)
!  INPUTS: odeint:
!  ystart= begin coordinates vector and speeds
!  nvar  = number of equations
!  t1    = begin integration time
!  t2    = end integration time
!          t2 may be > or < t1
!  eps   = absolute required precision for solution
!  h1    = time increment proposed at beginning (try value)
!  hmin  = minimum time increment
!
!  OUTPUTS:
!  ystart= end coordinates vector and speeds
!  h1    = final time increment
!  nok   = number of unchanged time steps
!  nbad  = number of modified time steps
!----------------------------------------------------------------------
!  DO WHILE (rsepmax .lt. rsep0 )
  DO WHILE (t1 .lt. tend_sec) 
     CALL odeint(Y,neq,t1,t2,tol,dt,dtf,nok,nbad)
     rsep     = SQRT((Y(5)-Y(1))**2+ (Y(7)-Y(3))**2)
     rsepmin  = MIN(rsep,rsepmin)
     rsepmax  = MAX(rsep,rsepmax)
     vs       = SQRT(Y(2)**2+Y(4)**2)
     engs     = 0.5D0*ms*vs**2
     engp     = 0.5D0*mp*(Y(6)**2+Y(8)**2)
     engt     = engs + engp
     angmtm    = ms*(y(1)*y(4) - y(3)*y(2)) + mp*(y(5)*y(8) - y(7)*y(6)) 
     rsrelx   = Y(1)-Y(5)    ! x location relative to planet
     rsrely   = Y(3)-Y(7)    ! y location relative to planet
     t1 = t2
     dt = dts
     t2 = t1 +10.*dt
     nsteps = nsteps+1
!     WRITE(*,11)t1,rsepmin,rsepmax,rsep
      WRITE(9,14)Y(1),Y(3),Y(5),Y(7),rsep,engs,vs,t1,rsrelx,rsrely,angmtm
!     WRITE(*,12)rsep,t1
!     WRITE(*,15)Y(1),Y(3),SQRT(Y(1)**2+Y(3)**2),rsep
!      WRITE(*,16)engs,Y(2),y(4),rsep,t1
!      WRITE(*,17)SQRT(y(1)**2+y(3)**2),t1
  ENDDO

  PRINT *,' '
  PRINT *,' At time = ', t2
  PRINT *,' external steps taken =',nsteps
  PRINT *,' no unchanged  steps =',nok
  PRINT *,' no modified time steps =',nbad
  DO i=1, neq
    WRITE(*,10) i, Y(i)
  END DO
  PRINT *,' '
  PRINT *,' Final time step = ', dtf
  STOP

 10 FORMAT('  Y(',I1,') = ',E13.6)
 11 FORMAT('time,rsepmin,rsepmax,rsep =',4(x,1pe12.4))
!12 FORMAT('xs,ys,xp,yp,rsep =',5(x,1pe12.4))
 14 FORMAT(11(x,1pe12.4))
 15 FORMAT('xs,ys,rs,rsep =',5(x,1pe12.4))
 16 FORMAT('engs,vsx,vys,rsep,time =',5(x,1pe12.4))
 17 FORMAT('rp,time =',5(x,1pe12.4))
END

!end of file drive.f90
