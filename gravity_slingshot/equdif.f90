!*****************************************************
!*              MODULE  Equdif.f90                   *
!* ------------------------------------------------- *
!*  Integration procedure odeint used by program     *
!*  tequdif.f90                                      *
!* ------------------------------------------------- *
!*             Author: Jean-Pierre Moreau, France    *
!*                      (www.jpmoreau.fr)            *
!*****************************************************
MODULE Equdif

CONTAINS

!=================================================================
Subroutine RK4(y, dydt, n, t, h, yout)
real*8 y(*), dydt(*), t, h, yout(*)
!=================================================================
!  From Pascal Version 1.1 - 12/31/1991 by J-P Dumont, France
!=================================================================
! This subroutine integrates the 1st order differential system:
!      dy/dt = F(y,t)
! by a fourth order Runge-Kutta method to advance the solution
! on interval h of the independant variable t:
!-----------------------------------------------------------------
! INPUTS:
! y     = State vector at begin of integration step
! dydt  = its derivative at the same point.
! n     = number of equations of system
! t     = i.s. value at begin of time step
! h     = time step
!-------------------------------------------------------------------
! OUTPUT:
! yout  = State vector at end of integration step.
!-------------------------------------------------------------------
! Programs using procedure RK4 must provide access to:
!
! SUBROUTINE Derivs(t,y,dydt)
! which returns the values of derivatives dydt at point t, knowing
! both t and the values of functions y.
!
!-------------------------------------------------------------------
      real*8 th, hh, h6
      real*8 dym(1:n), dyt(1:n), yt(1:n)

      hh = 0.5d0 * h
      h6 =h/6.d0
      th = t + hh
      do i = 1, n
	    yt(i) = y(i) + hh * dydt(i)
	  end do
      call Derivs(th,yt,dyt)
      do i = 1, n
	    yt(i) = y(i) + hh * dyt(i)
	  end do
      call Derivs(th,yt,dym)
      do i = 1, n
        yt(i) = y(i)+ h*dym(i)
        dym(i)= dyt(i) + dym(i)
      end do
      call Derivs(t+h,yt,dyt)
      do i = 1, n
	    yout(i) = y(i)+h6*(dydt(i)+dyt(i)+2.d0*dym(i))
	  end do
	  return
  End Subroutine RK4

!///////////////////////////////////////////////////////////////////////
Subroutine rkqc(y, dydt, n, t, htry, eps, yscal, hdid, hnext)
real*8 y(*), dydt(*), t, htry, eps, yscal(*), hdid, hnext
!=======================================================================
! From Pascal Version 1.2 - 03/26/1993 By J-P Dumont, France
!=======================================================================
! Runge-Kutta integration step with control of truncation local error 
! to obtain a required precision and adjust time step consequently
!-----------------------------------------------------------------------
! INPUTS:
! y      = State vector of size n
! dydt   = its derivative at begin value of independant variable, t
! n      = number of equations of system
! t      = begin value of independant variable
! htry   = time step proposed as a try
! eps    = precision requirement:
!          Max (ycalc[i] - yvrai[i])/yscal[i] < eps
! yscal  = normalization vector of solution.
!----------------------------------------------------------------------
! OUTPUTS:
!  y     = end state vector
!  t     = i.s. end value
!  hdid  = actual time step  
!  hnext = time step advised for the next integration step
!---------------------------------------------------------------------
! Programs using procedure RKQC must provide access to:
!
! Subroutine Derivs(t, y, dydt)
! which returns the values of derivatives dydt at point t, knowing
! both t and the values of functions y.
!
!=====================================================================
real*8 tsav,hh,h,temp,errmax
real*8 dysav(1:n), ysav(1:n), ytemp(1:n)
    
 	     pgrow=-0.20d0
         pshrnk=-0.25d0
         fcor=0.06666666d0   !1/15
         un = 1.d0;
         safety=0.9d0
         errcon=6.D-4
         tiny= 1.D-20    

         tsav= t          !Save begin time
         do i=1, n
           ysav(i) = y(i)
           dysav(i)= dydt(i)
         end do
         h = htry         !define increment for a try value
1        hh = 0.5d0*h     !take 2 half time steps
         call rk4(ysav,dysav,n,tsav,hh,ytemp)
         t = tsav + hh
         call Derivs(t,ytemp,dydt)
         call rk4(ytemp,dydt,n,t,hh,y)
         t = tsav + h
         if (t == tsav) then
           print *,' Pause in RKQC subroutine'
           print *,' Increment too small of independant variable'
           pause ' Press any key to continue...'
         end if
         call rk4(ysav,dysav,n,tsav,h,ytemp)
         errmax = 0.d0    !Evaluate error
         temp = 0.d0
         do i = 1, n
           ytemp(i) = y(i) - ytemp(i)    !ytemp = estimated error
           IF (yscal(i)>tiny)  temp = dabs(ytemp(i)/yscal(i))
           IF (errmax < temp)  errmax = temp
         end do
         errmax = errmax/eps             !real error / requirement
         if (errmax > un) then           !Error too big, reduce h
            h = safety*h*exp(pshrnk*dlog(errmax))
            GOTO 1                       !start again
         else                            !the step has been a success
            hdid = h                     !Calculate next time step
            IF (errmax > errcon) THEN
              hnext = safety*h*dexp(pgrow*dlog(errmax))
            ELSE 
			  hnext = 4.d0*h
            end if
         end if
         do i = 1, n 
		   y(i) = y(i) + ytemp(i)*fcor
         end do
		 return
End Subroutine rkqc

!//////////////////////////////////////////////////////////////////////
Subroutine odeint(ystart, nvar, t1, t2, eps, h1, hmin, nok, nbad)
real*8 ystart(*), t1, t2, eps, h1, hmin 
!================================================================
!  From Pascal Version 1.2 - 12/31/1991 By J-P Dumont, France
!================================================================
!  This subroutine integrates the 1st order differential system:
!       dy/dt = F(y,t)
!  where y and F are vectors of size nvar, between times
!  t1 and t2.
!
!  INPUTS:
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
!===================================================================
! Programs using procedure RKQC must provide access to:
!
! Subroutine Derivs(t, y, dydt)
! which returns the values of derivatives dydt at point t, knowing
! both t and the values of functions y.
!
!====================================================================
real*8 two,zero,tiny,tsav,t,hnext,hdid,h, yscal(1:nvar),y(1:nvar),dydt(1:nvar) 

     maxstp= 100000
     two   =    2
     zero  =    0
     tiny  = 1.d-20
	 dtsave = 0.d0       !not used here

     t= t1
     IF (t2 > t1) THEN 
	   h =   dabs(h1)
     ELSE 
	   h = - dabs(h1)
     END IF 
     nok = 0
     nbad = 0
     do i = 1, nvar 
	   y(i) = ystart(i)
     end do
     tsav = t - dtsave*two
     do nstp = 1, maxstp
         call Derivs(t,y,dydt)
         do i = 1, nvar 
		   yscal(i) = dabs(y(i))+dabs(dydt(i)*h)
         end do
         IF (((t+h-t2)*(t+h-t1)) > zero)   h = t2 - t
         call rkqc(y,dydt,nvar,t,h,eps,yscal,hdid,hnext)
         IF (hdid == h ) THEN 
		   nok = nok + 1
         ELSE 
		   nbad = nbad + 1
         END IF
         IF (((t-t2)*(t2-t1)) >= zero) THEN
           do i = 1, nvar 
		     ystart(i) = y(i)
           end do
           GOTO 99  !it is over
         END IF
         IF (abs(hnext) < hmin) THEN
           nok = -1   !error flag
           pause ' Time step too small!'
           GOTO 99
         END IF
         h = hnext
     end do
	 pause ' Pause in subroutine ODEINT - too many time steps!'
99   h1 = h
End Subroutine odeint

END MODULE EQUDIF

!end of file equdif.f90
